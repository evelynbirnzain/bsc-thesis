from tqdm import tqdm
import pygad
import numpy as np
import pandas as pd
import os
import io_utils


class GA:
    num_generations = 10
    num_parents_mating = 10
    sol_per_pop = 20
    init_range_low = 0
    parent_selection_type = "sss"
    keep_parents = 6
    crossover_type = "single_point"
    random_mutation_min_val = 0
    mutation_by_replacement = True
    mutation_type = "random"
    mutation_percent_genes = 20
    stop_criteria = ["saturate_30000"]
    instance = None

    def __init__(self, genes, scoring, plotting, init_pop_weights=None):
        self.genes = genes
        self.init_range_high = len(self.genes) - 1
        self.gene_space = range(0, len(self.genes))
        self.random_mutation_max_val = len(self.genes) - 1
        self.plotting = plotting
        self.scoring = scoring
        self.init_pop_weights = init_pop_weights

    def _fitness_func(self, solution, idx, key):
        selected_genes = self.genes[solution]
        return self.scoring.score(selected_genes, key)

    def _initial_population(self, num_genes):
        sols = []
        for _ in range(self.sol_per_pop):
            sel = list(self.init_pop_weights.sample(num_genes, weights=self.init_pop_weights).index)
            enc = np.isin(self.genes, sel).nonzero()[0]
            sols.append(enc)
        return sols

    def run_and_save_instance(self, fitness_func, num_genes, meth):
        with tqdm(total=self.num_generations) as pbar:
            self.instance = pygad.GA(
                num_generations=self.num_generations,
                num_parents_mating=self.num_parents_mating,
                fitness_func=fitness_func,
                sol_per_pop=self.sol_per_pop,
                num_genes=num_genes,
                init_range_low=self.init_range_low,
                init_range_high=self.init_range_high,
                parent_selection_type=self.parent_selection_type,
                keep_parents=self.keep_parents,
                mutation_type=self.mutation_type,
                mutation_percent_genes=self.mutation_percent_genes,
                gene_type=int,
                random_mutation_min_val=self.random_mutation_min_val,
                random_mutation_max_val=self.random_mutation_max_val,
                mutation_by_replacement=self.mutation_by_replacement,
                allow_duplicate_genes=False,
                on_generation=lambda _: pbar.update(1),
                save_best_solutions=True,
                stop_criteria=self.stop_criteria,
                suppress_warnings=True
            )
            if self.init_pop_weights:
                self.instance.initial_population = self._initial_population(num_genes)
            self.instance.run()

        self._save_results(meth)

    def _save_results(self, meth):
        best_sols = self.instance.best_solutions
        sols = pd.DataFrame(best_sols).apply(lambda x: self.genes[x])
        sols = io_utils.sort_genes(sols)
        sols = pd.concat([sols, sols.apply(self.scoring.score, axis=1)], axis=1).sort_values(by="single")

        run_dir = io_utils.create_run_dir(meth, self.instance.num_genes)
        sols.to_csv(os.path.join(run_dir, "best_sols.csv"), index=False)
        self.instance.plot_fitness().savefig(os.path.join(run_dir, "fit_vs_gen.svg"), format="svg")

        solution, solution_fitness, solution_idx = self.instance.best_solution()
        self.plotting.boxplot(self.genes[solution], os.path.join(run_dir, "best_sol.svg"))

        np.save(os.path.join(run_dir, "last_pop"), self.instance.population)

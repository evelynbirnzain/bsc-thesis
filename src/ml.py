from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from tqdm import tqdm
import pandas as pd
import os
from sklearn.model_selection import train_test_split
import random

import io_utils


class ML:
    def __init__(self, seed, X, y, scoring):
        random.seed(a=seed)
        self.feature_names = X.columns
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(X, y, test_size=0.2)
        self.scoring = scoring

    def _save_results(self, sols, run_dir, n_feat):
        states = sols["state"]
        sols = sols[range(n_feat)]
        sols = io_utils.sort_genes(sols)
        metrics = sols[range(n_feat)].apply(self.scoring.score, axis=1)
        sm = pd.concat([sols, metrics, states], axis=1).sort_values(by="single")
        sm.to_csv(os.path.join(run_dir, "best_sols.csv"), index=False)

    def fit_and_save(self, n_feat, n_rep):
        sols = []
        states = []
        run_dir = io_utils.create_run_dir("rf", n_feat)
        for i in tqdm(range(n_rep)):
            state = random.randint(1000, 100000)
            clf = RandomForestClassifier(n_jobs=5, random_state=state)
            sel = SelectFromModel(clf, max_features=n_feat)
            sel.fit(self.X_train, self.y_train)
            selected_feat = list(self.X_train.columns[(sel.get_support())])

            sols.append(selected_feat)
            states.append(str(state))

            if i % 10 == 0 or i == n_rep - 1:
                df_sols = pd.DataFrame(sols)
                df_sols["state"] = states
                self._save_results(df_sols, run_dir, n_feat)




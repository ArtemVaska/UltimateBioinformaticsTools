from concurrent.futures import ProcessPoolExecutor

import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier

SEED = 12345


class RandomForestClassifierCustom(BaseEstimator):
    def __init__(
        self, n_estimators=10, max_depth=None, max_features=None, random_state=SEED
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def generate_dt(self, args):
        tree_id, X, y = args

        n_samples, n_features = X.shape
        np.random.seed(self.random_state + tree_id)

        feature_idx = np.random.choice(
            n_features, size=self.max_features, replace=False
        )

        bootstrap_idx = np.random.choice(n_samples, size=n_samples, replace=True)
        X_bootstrap = X[bootstrap_idx]
        y_bootstrap = y[bootstrap_idx]

        tree = DecisionTreeClassifier(
            max_depth=self.max_depth,
            max_features=self.max_features,
            random_state=self.random_state + tree_id,
        )
        tree.fit(X_bootstrap[:, feature_idx], y_bootstrap)

        return tree, feature_idx

    def fit(self, X, y, n_jobs=1):
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            results = list(
                executor.map(
                    self.generate_dt,
                    [(tree_id, X, y) for tree_id in range(self.n_estimators)],
                )
            )

        self.trees = [results[i][0] for i in range(len(results))]
        self.feat_ids_by_tree = [results[i][1] for i in range(len(results))]

        return self

    def predict_proba(self, X, n_jobs=1):
        probas = 0
        for i, tree in enumerate(self.trees):
            probas += tree.predict_proba(X[:, self.feat_ids_by_tree[i]])
        avg_probas = probas / self.n_estimators

        return avg_probas

    def predict(self, X, n_jobs=1):
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)

        return predictions

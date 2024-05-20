from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier

SEED = 12345


class RandomForestClassifierCustom(BaseEstimator):
    """
    A custom implementation of a RandomForestClassifier with parallel realization.

    Attributes:
    [x] n_estimators: the number of decision trees in the forest
    [x] max_depth: the maximum depth of each decision tree
    [x] max_features: the number of features to consider when looking for the best split
    [x] random_state: the random seed for reproducibility

    Methods:
    [x] fit(X, y, n_jobs=n): fit the random forest on the input data using multiple processes n, default is 1
    [x] predict_proba(X, n_jobs=n): predict class probabilities for input data using multiple processes n, default is 1
    [x] predict(X, n_jobs=n): predict classes for input data using multiple processes n, default is 1
    """

    def __init__(
        self,
        n_estimators: int = 10,
        max_depth: int = None,
        max_features: int = None,
        random_state: int = SEED,
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def generate_dt(self, args: tuple):
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

    def fit(
        self, X: pd.Series, y: pd.Series, n_jobs: int = 1
    ) -> "RandomForestClassifierCustom":
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

    def loop_proba(self, args: tuple) -> pd.Series:
        X, i, tree = args
        proba = tree.predict_proba(X[:, self.feat_ids_by_tree[i]])
        return proba

    def predict_proba(self, X: pd.DataFrame, n_jobs: int = 1) -> float:
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            probas = executor.map(
                self.loop_proba, [(X, i, tree) for i, tree in enumerate(self.trees)]
            )

        avg_probas = sum(probas) / self.n_estimators
        return avg_probas

    def predict(self, X: pd.DataFrame, n_jobs=1) -> pd.Series:
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)

        return predictions

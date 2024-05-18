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

    def fit(self, X, y):
        self.classes_ = sorted(np.unique(y))
        n_samples, n_features = X.shape

        for i in range(self.n_estimators):
            np.random.seed(self.random_state + i)

            feature_idx = np.random.choice(n_features, size=self.max_features, replace=False)
            self.feat_ids_by_tree.append(feature_idx)

            bootstrap_idx = np.random.choice(n_samples, size=n_samples, replace=True)
            X_bootstrap = X[bootstrap_idx]
            y_bootstrap = y[bootstrap_idx]

            tree = DecisionTreeClassifier(
                max_depth=self.max_depth,
                max_features=self.max_features,
                random_state=self.random_state + i
            )
            tree.fit(X_bootstrap[:, self.feat_ids_by_tree[i]], y_bootstrap)
            self.trees.append(tree)
        return self

    def predict_proba(self, X):
        probas = 0
        for i, tree in enumerate(self.trees):
            probas += tree.predict_proba(X[:, self.feat_ids_by_tree[i]])
        avg_probas = probas / self.n_estimators
        return avg_probas

    def predict(self, X):
        probas = self.predict_proba(X)
        predictions = np.argmax(probas, axis=1)

        return predictions

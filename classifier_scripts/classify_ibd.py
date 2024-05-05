import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler
import numpy as np

def current_ibd(path):
    df = pd.read_csv(path, low_memory=False)

    # Select columns for microbiome data
    microbiome_data = df.columns[5:55]
    X = df[microbiome_data]
    y = df['Diagnosis']

    X_encoded = pd.get_dummies(X, columns=[col for col in X.columns if X[col].dtype == 'object'])

    X_encoded.fillna(X_encoded.mean(), inplace=True)

    X_encoded.replace([np.inf, -np.inf], np.nan, inplace=True)
    X_encoded.fillna(X_encoded.mean(), inplace=True)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_encoded)

    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

    classifier = RandomForestClassifier(random_state=42)
    classifier.fit(X_train, y_train)

    y_pred = classifier.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    print("###########FOR CURRENT PREDICTIONS#############")
    print(f"Accuracy: {accuracy:.2f}")

    scores = cross_val_score(classifier, X_scaled, y, cv=5)
    print(f"Average Cross-Validated Accuracy: {scores.mean():.2f}")

    feature_importances = classifier.feature_importances_
    importance_df = pd.DataFrame(feature_importances, index=X_encoded.columns, columns=['Importance'])
    importance_df = importance_df.sort_values(by='Importance', ascending=False)
    print("Top 20 features by importance:")
    for index, row in importance_df.head(60).iterrows():
        print(f"{index}: {row['Importance']:.6f}")

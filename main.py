import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler

def prepare_data(df):
    # Convert timestamp into a datetime object and sort data
    df['Timestamp'] = pd.to_datetime(df['Timestamp'])
    df.sort_values(by=['PatientID', 'Timestamp'], inplace=True)

    # Create lag features or differences that might indicate changes over time
    df['Previous_Diagnosis'] = df.groupby('PatientID')['Diagnosis'].shift(1)
    df.fillna({'Previous_Diagnosis': 'No Diagnosis'}, inplace=True)  # Fill NA for first entry of each patient

    return df

def current_ibd():
    # Load data
    df = pd.read_csv('LiuLloydGevers1_AbudanceData_WSpecies_WMetadata_IBD_Normalized_WStudy.csv', low_memory=False)

    # Select columns for microbiome data
    microbiome_data = df.columns[5:55]
    X = df[microbiome_data]
    y = df['Diagnosis']

    # Identify and one-hot encode categorical columns
    X_encoded = pd.get_dummies(X, columns=[col for col in X.columns if X[col].dtype == 'object'])

    # Handle missing values
    X_encoded.fillna(X_encoded.mean(), inplace=True)

    # Check and handle infinite values
    X_encoded.replace([np.inf, -np.inf], np.nan, inplace=True)
    X_encoded.fillna(X_encoded.mean(), inplace=True)

    # Scale data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_encoded)

    # Split data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)

    # Train classifier
    classifier = RandomForestClassifier(random_state=42)
    classifier.fit(X_train, y_train)

    # Evaluate classifier
    y_pred = classifier.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    print(f"Accuracy: {accuracy:.2f}")

    # Cross-validation
    scores = cross_val_score(classifier, X_scaled, y, cv=5)
    print(f"Average Cross-Validated Accuracy: {scores.mean():.2f}")

    # Get feature importances
    feature_importances = classifier.feature_importances_
    importance_df = pd.DataFrame(feature_importances, index=X_encoded.columns, columns=['Importance'])
    importance_df = importance_df.sort_values(by='Importance', ascending=False)
    print("Top 20 features by importance:")
    for index, row in importance_df.head(21).iterrows():
        print(f"{index}: {row['Importance']:.6f}")

if __name__ == '__main__':
    current_ibd()
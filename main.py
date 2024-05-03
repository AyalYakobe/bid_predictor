from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.stattools import adfuller


def predict_microbiome_state(filepath, target_column, date_column, order=(1, 1, 1), forecast_steps=5):

    # Load data with specified dtype for problematic columns if known
    # Example: df = pd.read_csv(filepath, dtype={'Column1': float, 'Column6': float})
    df = pd.read_csv(filepath, low_memory=False)  # Set low_memory to False if unsure

    # Convert date column to datetime and ensure it's sorted
    df[date_column] = pd.to_datetime(df[date_column])
    df.sort_values(by=date_column, inplace=True)  # Ensure data is sorted
    df.set_index(date_column, inplace=True)

    # Setting a frequency for the datetime index (e.g., 'W' for weekly data)
    df.index = pd.DatetimeIndex(df.index).to_period('W')  # Change 'W' based on your data frequency

    # Select the target variable
    target = df[target_column].dropna()

    # Check for stationarity
    result = adfuller(target)
    print('ADF Statistic: %f' % result[0])
    print('p-value: %f' % result[1])

    # Differencing if necessary
    if result[1] > 0.05:
        target = target.diff().dropna()

    # Fit ARIMA model
    model = ARIMA(target, order=order)
    model_fit = model.fit()

    # Forecast
    forecast = model_fit.forecast(steps=forecast_steps)

    return forecast


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
    print("###########FOR CURRENT PREDICTIONS#############")
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
    #current_ibd()
    # future_microbiome_state()

    file_path = 'LiuLloydGevers1_AbudanceData_WSpecies_WMetadata_IBD_Normalized_WStudy.csv'
    target_column = 'Bacteria;Firmicutes;Clostridia'
    date_column = 'Collection Week'
    forecast = predict_microbiome_state(file_path, target_column, date_column)
    print(forecast)
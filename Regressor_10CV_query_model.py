# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 14:56:16 2024

@author: aleniak
"""

import pandas as pd
import joblib
import sys
import os
import math
import numpy as np
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

# Function to make predictions and evaluate the model
def predict_and_evaluate(csv_file, model_file):
    try:
        # Load the data from the CSV file
        data = pd.read_csv(csv_file)  # Load the data where the first column is the actual values, and the rest are features
        X = data.iloc[:, 1:]  # Select all columns except the first one as input data (features)
        y_actual = data.iloc[:, 0]  # The first column is the actual values (target)

        # Load the trained model
        model = joblib.load(model_file)  # Load the trained model using joblib

        # Make predictions based on the input features
        y_pred = model.predict(X)

        # List to store the results for each row
        results = []

        # Lists to store global metrics (RMSE and MAE for each row)
        rmse_values = []
        mae_values = []

        # Iterate through actual and predicted values
        for actual, predicted in zip(y_actual, y_pred):
            # Calculate metrics for each row
            mse = mean_squared_error([actual], [predicted])
            rmse = math.sqrt(mse)  # Calculate RMSE
            mae = mean_absolute_error([actual], [predicted])  # Calculate MAE

            # Round values to two decimal places
            predicted_rounded = round(predicted, 2)
            rmse_rounded = round(rmse, 3)
            mae_rounded = round(mae, 3)

            # Store the results for each row
            results.append({
                'Actual Value': round(actual, 3),  # Round the actual value
                'Predicted Value': predicted_rounded,  # Round the predicted value
                'RMSE': rmse_rounded,  # RMSE for a single row
                'MAE': mae_rounded  # MAE for a single row
            })

            # Collect RMSE and MAE for global metrics
            rmse_values.append(rmse)
            mae_values.append(mae)

        # Calculate global metrics: average and standard deviation for RMSE and MAE
        averaged_rmse = round(np.mean(rmse_values), 3)
        std_dev_rmse = round(np.std(rmse_values), 3)
        averaged_mae = round(np.mean(mae_values), 3)
        std_dev_mae = round(np.std(mae_values), 3)
        r2 = round(r2_score(y_actual, y_pred), 4)  # Calculate R² score

        # Add global metrics at the end of the results
        results.append({
            'Actual Value': 'Global Metrics',
            'Predicted Value': '',
            'RMSE': '',
            'MAE': ''
        })
        results.append({
            'Actual Value': '',
            'Predicted Value': '',
            'RMSE': f'Average RMSE: {averaged_rmse}',
            'MAE': f'Average MAE: {averaged_mae}'
        })
        results.append({
            'Actual Value': '',
            'Predicted Value': '',
            'RMSE': f'Std Dev RMSE: {std_dev_rmse}',
            'MAE': f'Std Dev MAE: {std_dev_mae}'
        })
        results.append({
            'Actual Value': '',
            'Predicted Value': '',
            'RMSE': f'R2: {r2}',
            'MAE': ''
        })

        # Save the results to a CSV file
        result_df = pd.DataFrame(results)  # Convert the results to a DataFrame
        output_file = os.path.splitext(csv_file)[0] + "_predictions.csv"  # Create the output file name
        result_df.to_csv(output_file, index=False)  # Save the results to a file
        print(f"Predictions saved to {output_file}")  # Success message for file save

    except Exception as e:
        print(f"Error: {e}")  # Error handling

# Main script logic
if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print(f"Usage: python {os.path.basename(__file__)} <csv_file> <model_file>")
        sys.exit(1)

    # Get the CSV file and model file paths from command-line arguments
    csv_file = sys.argv[1]
    model_file = sys.argv[2]

    # Call the function to make predictions and evaluate the model
    predict_and_evaluate(csv_file, model_file)

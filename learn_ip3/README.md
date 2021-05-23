# Train models generalizing in IP3

This directory trains a model in the IP3 concentration for a single volume and number of IP3Rs.

It uses the data in the [ml_training_data](../stochastic_simulations/ml_training_data) directory - visit that before this.

Run it in order following this guide.

## Params from data

Either:
* Pull the pre-loaded data using `dvc pull`, or
* Run the notebook `cache_params0.nb` to identify the ML parameters from the data, and write them to the `cache` folder.

## Derivatives

Either:
* Pull the pre-loaded data using `dvc pull`, or
* Run the notebook `cache_derivs0.nb` to differentiate the ML parameters using total variation regularization, and write them to the `cache` folder.

## Make the training data

Head to the [train/training_data](train/training_data) folder.
1. Create the graphs with the `make_graph_rxns.nb` notebook for the reaction-based model. To run this, **first** run the initialization cells in the `funcs_network_rxns.nb` to make the necessary function definitions.
2. Similarly, you can create the graphs for the parameter-based model with the corresponding methods.
3. Create the training data with the `make_training_data.nb` notebook.

The results are in the `training_data_params` and `training_data_rxns` folders.

## Train the model

Head to the [train/train](train/train) folder.

Either:
* Pull the pre-loaded trained networks using `dvc pull`, or
* Use the notebook `train.nb` to train the models.

The resulting trained networks and data are in the `trained` folder.

## Extra: transformations Jacobian (supplemental material)

To calculate the Jacobians for the parameter transformations, see the [transformation_jacobian](transformation_jacobian) folder.
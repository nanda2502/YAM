# YAM: Yet Another Model

## Project Description

**YAM (Yet Another Model)** is a Haskell-based computational tool designed to study the dynamics of trait learning. The model calculates the expected number of steps to absorption (i.e., when all traits are learned) under different learning strategies. 

Given a number of nodes, YAM generates all possible adjacency matrices representing learning pathways. For each adjacency matrix, the model applies both Random Learning and Payoff-Based Learning strategies across different alpha parameters, calculating the associated transition matrices and expected steps to absorption. The results are saved as CSV files for further analysis.

## Installation and Running

### Prerequisites

- GHC (The Glasgow Haskell Compiler)
- Cabal (Common Architecture for Building Applications and Libraries)

### Steps

1. **Clone the Repository:**

    ```bash
    git clone <repository-url>
    cd yam
    ```

2. **Build the Project:**

    ```bash
    cabal build
    ```

3. **Run the Model:**

    ```bash
    cabal run YAM [number_of_nodes] [save_transition_matrices]
    ```

    - `number_of_nodes`: The desired number of nodes (default is 4 if not specified)
    - `save_transition_matrices`: A boolean flag (True/False) to control whether transition matrices are saved (default is True if not specified)

    Examples:
    - To run the model with 4 nodes and save transition matrices (default behavior):
      ```bash
      cabal run YAM 4
      ```
    - To run the model with 6 nodes and not save transition matrices:
      ```bash
      cabal run YAM 6 False
      ```
    - To run the model with default 4 nodes and not save transition matrices:
      ```bash
      cabal run YAM 4 False
      ```

The program will generate an `output` directory, saving the summary of expected steps to absorption in a CSV file. If `save_transition_matrices` is set to True, it will also save individual transition matrices as CSV files. Note that the `output` directory is cleared at the start of each run.

## Output

- `expected_steps.csv`: A summary of expected steps to absorption for each combination of adjacency matrix, learning strategy, and alpha value.
- Transition matrix files (if saving is enabled): Named as `transition_mat_<adjacency_matrix>_strategy_<strategy>_alpha_<alpha>.csv`

## Note

Saving transition matrices can significantly increase the execution time and storage requirements, especially for larger numbers of nodes. Use the `save_transition_matrices` flag to control this behavior based on your needs.

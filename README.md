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
    cabal run YAM <number_of_nodes>
    ```

    Replace `<number_of_nodes>` with the desired number of nodes. For example, to run the model with 4 nodes:

    ```bash
    cabal run YAM 4
    ```

The program will generate an `output` directory, saving transition matrices and summarizing expected steps to absorption in CSV files.

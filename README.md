## Model Overview

This codebase implements a Markov chain model of social learning in networks of interdependent traits. Individuals learn traits from others in a population, where traits have prerequisite dependencies (e.g., you must learn basic arithmetic before calculus). The model compares how different social learning strategies perform under varying levels of constraints on learning.

### How It Works

1. **Trait Networks**: Represented as directed acyclic graphs where nodes are traits and edges define prerequisites
2. **Learning Strategies**: Individuals use different strategies to select which traits to attempt learning (payoff-biased, conformist, proximal, etc.)
3. **Markov Analysis**: The system evolves as a Markov chain over repertoire states (which traits each learner knows)
4. **Metrics**: We track performance (cumulative payoff), learning efficiency (successful transitions), and population diversity (variation in repertoires)

### Program Flow
```
main.cpp 
  ↓
ModelConfig::makeCombinations()  // Generate parameter sweeps
  ↓
ModelConfig::processRepl()       // Run one parameter combination
  ↓
MarkovChain::computeMarkovChain() // Build transition matrix, compute metrics
  ↓
IO::writeAndCompressCSV()        // Save results
```

**Input**: Adjacency matrices in `data/adj_mat_<postfix>.csv` (one per line, flattened)  
**Output**: Results in `output/yam_out_<adj_idx>.csv` with columns for each metric and parameter

---

## Core Modules

### ModelConfig (ModelConfig.cpp/hpp)
**Purpose**: Orchestrates parameter combinations and manages simulation runs

**Key Functions**:
- `makeCombinations()`: Generates all parameter combinations tested in the paper
- `processRepl()`: Executes one replication across payoff permutations
- `makeShuffles()`: Generates all (or sampled) permutations of payoff assignments
- `computeTransitiveClosure()`: Adds edges between all ancestor-descendant pairs

**Parameter Combinations by Figure**:

**Base cases** (default: slope=2, transparency=0, edge_weight=1.0, alpha=0, payoffdist=0, distribution=Learnability):
- **Main text**: Figure 1, Figure 2A-C
- **Supplementary**: Figure S2 (time series), S3 (distributions), S6 (learning time), S7A-C (absolute performance), S11 (Markov vs simulation), Table S1

**Sensitivity tests** (n=8 only):
- **Figure 2D-E**: Trait expression bias (`distribution = Depth/Payoffs`)
- **Figure 2F**: Payoff-depth correlation (`alpha = 1.0`)
- **Figure 2G-I**: Soft dependencies (`edge_weight < 1.0`, `closure = 1`)
- **Figure 2J-L**: Transparency (`transparency = 1, 10, 100`)
- **Figure S4**: Varying slopes (`slope = 1.25, 2.0, 5.0`)
- **Figure S5**: Prestige2 implementation (repertoire size instead of total payoff)
- **Figure S7D-L**: Absolute performance versions of Figure 2 sensitivity tests
- **Figure S8**: Skewed payoffs (`payoffdist = 1`)
- **Figure S10**: Smaller networks (requires separate runs with `postfix = 3, 4, 5`)
- **Figure S12**: Transitive reduction vs closure (`closure = 0` vs `closure = 1` with varying `edge_weight`)

**Parameters you can vary**:
- `slope`: Strategy bias strength (default: 2.0)
- `transparency`: Knowledge of prerequisite structure (0 = opaque, higher = more transparent)
- `edgeWeight`: Softness of dependencies (1.0 = deterministic, <1.0 = probabilistic)
- `distribution`: How traits are expressed (Learnability, Depth, Payoffs, etc.)
- `alpha`: Correlation between payoffs and trait depth (0 = random, 1 = correlated)
- `payoffDist`: Payoff distribution (0 = uniform spacing, 1 = one high + rest low)
- `closure`: Use transitive closure (0 = minimal edges, 1 = all ancestor edges)

---

### MarkovChain (MarkovChain.cpp/hpp)
**Purpose**: Constructs and analyzes the Markov chain over repertoire states

**Key Functions**:
- `computeMarkovChain()`: Main analysis pipeline
  1. Generate initial trait/state frequencies (first pass with random learner)
  2. Bias frequencies according to `TraitDistribution` parameter
  3. Build final transition matrix for the specified strategy
  4. Compute all metrics (payoff, transitions, variation, absorption time)

- `buildTransitionMatrix()`: Assembles state-to-state transition probabilities
- `biasTraitFrequencies()`: Adjusts observability of traits (for Fig 2D-E)
- `reorderTransitionMatrix()`: Separates transient from absorbing states
- `computeExpectedPayoffAtNSteps()`: Evolves state distribution over 20 steps
- `computeExpectedTimeToAbsorption()`: Uses fundamental matrix to find absorption time
- `computeStationaryVariation()`: Measures diversity in quasi-stationary distribution

### Metrics Tracked Over 20 Steps
- **Payoff**: Expected cumulative reward
- **Transitions**: Expected number of successful learning attempts
- **Variation**: Average Jaccard distance between repertoires in the population
- **Absorption time**: Expected steps until all traits are learned (computed once)

---

### Learning (Learning.cpp/hpp)
**Purpose**: Defines trait selection logic for each learning strategy

**Key Functions**:
- `learnability()`: Computes probability each trait can be learned from current repertoire
  - Deterministic: 1.0 if all prerequisites known, 0.0 otherwise
  - Probabilistic: Product of (1 - edgeWeight) for each missing prerequisite

- `normalizedWeights()`: Combines strategy bias with learnability
  - Filters out already-known traits and unobservable traits
  - Applies perceived learnability penalty based on `transparency`
  - Returns probability distribution over traits to attempt

- `baseWeights()`: Strategy-specific weighting before normalization
  - **Payoff**: `frequency × payoff^slope`
  - **Conformity**: `frequency^slope`
  - **Prestige**: `Σ(demonstrator_frequency × demonstrator_payoff^slope)` for each trait
  - **Proximal**: `Σ(demonstrator_frequency × Δ^(-slope))` where Δ = repertoire gap
  
- `generateReachableRepertoires()`: BFS through state space from naive state
- `transitionFromState()`: Computes transitions from one repertoire to all reachable neighbors

---

## Supporting Modules

### Graph (Graph.cpp/hpp)
- `computeDistances()`: BFS from root to find depth of each trait (for payoff correlation tests)
- `parentTraits()`: Returns direct prerequisites of a trait

### Payoffs (Payoffs.cpp/hpp)
- `generatePayoffs()`: Creates payoff vector for non-root traits
  - `alpha=0`: Payoffs permuted randomly (most analyses)
  - `alpha=1`: Payoffs increase with trait depth (Figure 2F)
  - `payoffDist>0`: High/low skew distributions (Figure S8)

### LinAlg (LinAlg.cpp/hpp)
- `decomposeLU()`: LU factorization with partial pivoting
- `solveUsingLU()`: Solve linear systems (used for fundamental matrix)

### IO (IO.cpp/hpp)
- `readAdjacencyMatrices()`: Parses three formats:
  1. Binary: `"01001100..."` → unweighted graph
  2. Compact weighted: `"03506920..."` → each digit is 0.0–0.9
  3. Comma-separated: `"0.0,0.3,0.5,..."` → full precision weights
- `writeAndCompressCSV()`: Outputs results
- `adjMatrixToFlattenedString()`: Converts matrix back to string for storage

### StringUtils & Types
- **StringUtils**: Formatting utilities (strategy names, CSV rows, state strings)
- **Types.hpp**: Core type definitions
  - `Strategy`: {Random, Payoff, Proximal, Prestige, Conformity, Prestige2}
  - `TraitDistribution`: {Learnability, Uniform, Depth, Payoffs}
  - `ParamCombination`: Struct bundling all parameters for one run

---

## Running the Model

### Basic Usage
```bash
./yam <adj_matrix_index> [<postfix>]
```

**Example**: Run analysis on the 42nd adjacency matrix from the 8-trait dataset
```bash
./yam 42 8
```

This will:
1. Load adjacency matrices from `data/adj_mat_8.csv`
2. Process the matrix at index 42 with all parameter combinations defined in `makeCombinations()`
3. Write results to `output/yam_out_42.csv`

### Understanding the Output

Each row in the output CSV represents one time step of one parameter combination:
```csv
num_nodes,adj_mat,alpha,strategy,repl,steps,step_payoff,step_transitions,step_variation,slope,distribution,absorbing,payoffdist,edge_weight,transparency,closure
8,01000000...,0.0,Proximal,0,5,3.42,4.1,0.23,Learnability,12.5,0.18,0,1.0,0.0,0
```

**Key columns**:
- `adj_mat`: Flattened adjacency matrix (64 chars for 8×8)
- `strategy`: Which learning strategy
- `steps`: Time step (1-20)
- `step_payoff`: Expected cumulative payoff at this step
- `slope`: Bias strength parameter
- `transparency`: Structural knowledge parameter
- `edge_weight`: Softness of prerequisites
- `absorbing`: Expected time to learn all traits (reported at step 1)



---

## Common Modifications

### Adding a New Learning Strategy

1. Add enum to `Types.hpp`: `enum Strategy { ..., MyNewStrategy };`
2. Implement weighting in `Learning.cpp`: Add case to `baseWeights()`
3. Add slope defaults in `ModelConfig.cpp`: `returnSlopeVector()`
4. Add string conversion in `StringUtils.cpp`: `strategyToString()`

### Testing a New Network Structure

1. Create adjacency matrix file: `data/adj_mat_custom.csv`
   - One flattened matrix per line (binary or weighted format)
2. Run with your postfix: `./yam 0 custom`
3. If desired, add special case to `generatePayoffs()` in `Payoffs.cpp` to set custom payoffs for your graph.

### Varying a New Parameter

1. Add field to `ParamCombination` struct in `Types.hpp`
2. Update `makeCombinations()` in `ModelConfig.cpp` to generate combinations
3. Pass parameter through `computeMarkovChain()` and relevant functions
4. Add output column in `formatResults()` in `StringUtils.cpp`

---

## Build Instructions

### Windows

To build and run YAM on Windows:

1. **Install MSYS2**

   Download and install MSYS2 from [https://www.msys2.org/](https://www.msys2.org/).

2. **Install Build Tools**

   Open the **MSYS2 MinGW 64-bit** terminal (not the MSYS terminal) and update the package database:

   ```bash
   pacman -Syu
   ```

   Close and reopen the terminal as prompted.

   Install the required packages:

   ```bash
   pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-cmake mingw-w64-x86_64-make
   ```

3. **Clone the Repository**

   Navigate to your preferred directory and clone the repository:

   ```bash
   git clone https://github.com/nanda2502/YAM.git
   ```

4. **Build the Project**

   ```bash
   cd YAM
   mkdir build
   cd build
   cmake -G "MinGW Makefiles" ..
   mingw32-make
   ```

5. **Run the Program**

   The executable `yam.exe` will be generated in the `build` directory. Run it with:

   ```bash
   ./yam <number_of_nodes>
   ```

   For example, to run the program with 5 nodes:

   ```bash
   ./yam 5
   ```

### Linux and macOS

To build and run YAM on Linux or macOS:

1. **Install Build Tools**

   Ensure you have the following installed:

   - **C++ Compiler**: GCC or Clang
   - **CMake**
   - **Make**
   - **Git**

   On **Ubuntu/Debian**, you can install them with:

   ```bash
   sudo apt update
   sudo apt install build-essential cmake make git
   ```

   On **macOS**, using Homebrew:

   ```bash
   brew install cmake git
   ```

2. **Clone the Repository**

   ```bash
   git clone https://github.com/nanda2502/YAM.git
   ```

3. **Build the Project**

   ```bash
   cd YAM
   mkdir build
   cd build
   cmake ..
   make
   ```

4. **Run the Program**

   The executable `yam` will be generated in the `build` directory. Run it with:

   ```bash
   ./yam <number_of_nodes>
   ```

   For example, to run the program with 5 nodes:

   ```bash
   ./yam 5
   ```
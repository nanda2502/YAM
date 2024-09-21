{-# LANGUAGE FlexibleContexts #-}
module Main where

import Types
import Graph
import Payoff (generatePayoffs)
import Learning
import LinearAlgebra (solveLinearSystem)
import Utils (writeMatrixToCSV, formatResult, readAdjacencyMatrices)

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector as V
import qualified Data.Map.Strict as M
import System.Random (StdGen, getStdGen)
import Control.Monad (forM, when)
import System.Environment (getArgs)
import System.Directory (doesDirectoryExist, removeDirectoryRecursive, createDirectory)
import System.FilePath ((</>))
import Text.Printf (printf)

-- | Computes the expected number of steps to absorption for a given strategy and alpha
computeExpectedSteps :: AdjacencyMatrix -> Strategy -> Double -> StdGen -> IO (Double, [[Double]])
computeExpectedSteps adjacencyMatrix strategy alpha gen = do
    let rootNode = 0  -- Root node index
    let distances = computeDistances adjacencyMatrix rootNode

    -- Generate payoffs based on distances and alpha
    let (payoffs, _) = generatePayoffs distances alpha gen

    -- Number of traits
    let n = V.length adjacencyMatrix

    -- Generate all reachable repertoires
    let repertoiresList = generateReachableRepertoires strategy adjacencyMatrix payoffs

    -- Assign indices to repertoires
    let repertoiresWithIndices = zip [0..] repertoiresList

    -- Create a mapping from repertoire to index
    let repertoireIndexMap = M.fromList [ (VU.toList r, i) | (i, r) <- repertoiresWithIndices ]

    -- Total number of states
    let numStates = length repertoiresList

    -- Build the transition matrix
    let transitionMatrix = V.generate numStates $ \i ->
            let repertoire = repertoiresList !! i
                transitions = transitionFromState strategy repertoire adjacencyMatrix payoffs
                stayProb = stayProbability strategy repertoire adjacencyMatrix payoffs
                -- Map repertoires to indices and probabilities
                transitionIndices = [ (repertoireIndexMap M.! VU.toList r', prob) | (r', prob) <- transitions ]
                -- Collect all transitions, including stay probability
                transitionsWithStay = (i, stayProb) : transitionIndices
                -- Sum up probabilities per index using Map
                probMap = M.fromListWith (+) transitionsWithStay
                -- Create the row vector by mapping the probabilities from probMap
                rowList = [ M.findWithDefault 0.0 idx probMap | idx <- [0..numStates -1] ]
                rowWithTransitions = V.fromList rowList
            in rowWithTransitions

    -- Converting transitionMatrix to list of lists
    let transitionMatrixList = map V.toList $ V.toList transitionMatrix

    -- Identify absorbing and transient states
    let isAbsorbingState = VU.all id  -- All traits are learned
        absorbingStates = [ (i, r) | (i, r) <- repertoiresWithIndices, isAbsorbingState r ]
        transientStates = [ (i, r) | (i, r) <- repertoiresWithIndices, not (isAbsorbingState r) ]

    -- Reorder the state indices: transient states first, then absorbing states
    let reorderedStateIndices = map fst transientStates ++ map fst absorbingStates

    -- Create a mapping from old indices to new indices
    let oldToNewIndexMap = M.fromList $ zip reorderedStateIndices [0..]

    -- Reorder the transition matrix rows
    let reorderedTransitionMatrixRows = [ transitionMatrixList !! oldIndex | oldIndex <- reorderedStateIndices ]

    -- Reorder the columns of each row
    let reorderedTransitionMatrix = [ [ row !! oldIndex | oldIndex <- reorderedStateIndices ] | row <- reorderedTransitionMatrixRows ]

    -- Number of transient and absorbing states
    let numTransientStates = length transientStates
        _numAbsorbingStates = length absorbingStates

    -- Extract Q and R matrices
    let qMatrix = [ take numTransientStates row | row <- take numTransientStates reorderedTransitionMatrix ]
        -- rMatrix = [ drop numTransientStates row | row <- take numTransientStates reorderedTransitionMatrix ]

    -- Build (I - Q)
    let identityMatrix = [ [ if i == j then 1.0 else 0.0 | j <- [0..numTransientStates -1]] | i <- [0..numTransientStates -1] ]
        iMinusQ = zipWith (zipWith (-)) identityMatrix qMatrix

    -- Build vector b (ones vector)
    let bVector = replicate numTransientStates 1.0

    -- Solve (I - Q) * t = 1 using solveLinearSystem
    let maybeTSolution = solveLinearSystem iMinusQ bVector

    tSolution <- case maybeTSolution of
                    Nothing -> do
                        -- Save the matrix to a CSV file
                        writeMatrixToCSV "matrix_failure.csv" iMinusQ
                        -- Then output the error
                        error "Failed to solve linear system: Matrix is singular or not invertible"
                    Just t -> return t

    -- Find the index of the initial state (root trait learned)
    let initialRepertoire = VU.accum (const id) (VU.replicate n False) [(0, True)]  -- Root trait is learned
        initialStateIndex = repertoireIndexMap M.! VU.toList initialRepertoire
        initialStateNewIndex = oldToNewIndexMap M.! initialStateIndex

    -- Expected number of steps from the initial state
    let expectedSteps = tSolution !! initialStateNewIndex

    return (expectedSteps, transitionMatrixList)

main :: IO ()
main = do
    args <- getArgs
    let (numNodes, saveTransitionMatrices) = parseArgs args
        n = numNodes

    let alphas = [0.0, 1.0, 2.0]  -- Alpha parameters
        strategies = [RandomLearning, PayoffBasedLearning]  -- Learning strategies

    gen <- getStdGen  -- Use system's standard random generator

    let outputDir = "./output"

    dirExists <- doesDirectoryExist outputDir 
    when dirExists $ removeDirectoryRecursive outputDir
    createDirectory outputDir

    -- Generate all adjacency matrices
    adjacencyMatrices <- readAdjacencyMatrices n

    putStrLn $ "Loaded " ++ show (length adjacencyMatrices) ++ " unique adjacency matrices." 

    -- For each adjacency matrix, compute expected steps for each alpha and strategy
    results <- forM adjacencyMatrices $ \adjMatrix -> do
        let adjMatrixBinary = adjMatrixToBinary adjMatrix
        res <- forM strategies $ \strategy -> do
            forM alphas $ \alpha -> do
                (expectedSteps, transitionMatrix) <- computeExpectedSteps adjMatrix strategy alpha gen

                let strategyStr = case strategy of
                                    RandomLearning      -> "RandomLearning"
                                    PayoffBasedLearning -> "PayoffBasedLearning"

                -- Save the transition matrix if the flag is set
                when saveTransitionMatrices $ do
                    let alphaStr = printf "%.2f" alpha  -- Format alpha with two decimal places
                        fileName = "transition_mat_" ++ adjMatrixBinary ++ "_strategy_" ++ strategyStr ++ "_alpha_" ++ alphaStr ++ ".csv"
                        filePath = outputDir </> fileName
                    writeMatrixToCSV filePath transitionMatrix

                return (n, adjMatrixBinary, alpha, strategyStr, expectedSteps)
        return (concat res)

    -- Flatten the results
    let flatResults = concat results

    -- Prepare CSV data with header
    let csvHeader = "NumNodes,AdjacencyMatrix,Alpha,Strategy,ExpectedSteps"
        csvData = csvHeader : map formatResult flatResults

    -- Write results to CSV
    let outputCsvPath = outputDir </> "expected_steps.csv"
    writeFile outputCsvPath $ unlines csvData

    putStrLn $ "Expected steps to absorption saved to '" ++ outputCsvPath ++ "'"

parseArgs :: [String] -> (Int, Bool)
parseArgs args = case args of
    [] -> (4, False)  -- Default values
    [n] -> (read n, True)
    [n, arg2] 
        | isBoolString arg2 -> (read n, read arg2)
        | otherwise         -> error "Second argument must be True or False."
    _ -> error "Usage: program [numNodes] [saveTransitionMatrices]"
  where
    isBoolString s = s == "True" || s == "False"
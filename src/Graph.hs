{-# LANGUAGE FlexibleContexts #-}
module Graph (
    computeDistances,
    parentTraits,
    generateAdjacencyMatrices,
    adjMatrixToBinary
) where

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Unboxed.Mutable as VUM
import qualified Data.Vector as V
import Types (Trait, AdjacencyMatrix)


-- | Computes distances from the root node (0) to all other nodes
computeDistances :: AdjacencyMatrix -> Trait -> VU.Vector Int
computeDistances adjMatrix root = VU.create $ do
    let n = V.length adjMatrix
    distances <- VUM.replicate n (-1)  -- Initialize distances to -1
    let loop [] = return ()
        loop ((node, dist):rest) = do
            d <- VUM.read distances node
            if d == -1
                then do
                    VUM.write distances node dist
                    let neighbors = [ i | i <- [0 .. n - 1], (adjMatrix V.! node) VU.! i ]
                    loop (rest ++ [ (neighbor, dist + 1) | neighbor <- neighbors ])
                else
                    loop rest
    loop [ (root, 0) ]
    return distances

-- | Gets parent traits of a given trait
parentTraits :: AdjacencyMatrix -> Trait -> [Trait]
parentTraits adjMatrix trait =
    [ i
    | i <- [0 .. V.length adjMatrix - 1]
    , (adjMatrix V.! i) VU.! trait
    ]

-- | Generates all possible adjacency matrices under the given rules
generateAdjacencyMatrices :: Int -> IO [AdjacencyMatrix]
generateAdjacencyMatrices n = do
    let allAdjMatrices = [ buildAdjacencyMatrix n parentChoices
                            | parentChoices <- parentChoicesCombos,
                              node0HasChildren parentChoices ]
    return allAdjMatrices
  where
    parentChoicesCombos = sequence [ [0..i-1] | i <- [1..n-1] ]

    -- Build adjacency matrix from parentChoices
    buildAdjacencyMatrix :: Int -> [Int] -> AdjacencyMatrix
    buildAdjacencyMatrix n parentChoices = adjacencyMatrix where
        edges = [ (p, i) | (i, p) <- zip [1..(n-1)] parentChoices ]

        hasEdge i j = (i, j) `elem` edges

        adjacencyMatrix = V.fromList [ VU.fromList [ hasEdge i j | j <- [0 .. n - 1] ] | i <- [0 .. n - 1] ]

    node0HasChildren parentChoices = 0 `elem` parentChoices

-- | Convert adjacency matrix to a binary string
adjMatrixToBinary :: AdjacencyMatrix -> String
adjMatrixToBinary adjMatrix =
    let boolList = concatMap VU.toList (V.toList adjMatrix)
    in map (\b -> if b then '1' else '0') boolList
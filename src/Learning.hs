{-# LANGUAGE FlexibleContexts #-}
module Learning (
    learnability,
    baseWeights,
    weightVector,
    transitionFromState,
    stayProbability,
    generateReachableRepertoires
) where

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector as V
import Types (Strategy(..), Repertoire, AdjacencyMatrix, PayoffVector)
import Graph (parentTraits)

-- | Computes learnability of traits given the current repertoire and adjacency matrix
learnability :: Repertoire -> AdjacencyMatrix -> VU.Vector Bool
learnability repertoire adjMatrix =
    VU.generate (VU.length repertoire) $ \j ->
        not (repertoire VU.! j) && -- Not already learned
        all (repertoire VU.!) (parentTraits adjMatrix j) -- All parents are learned

-- | Base weights assigned by the strategy to each trait
baseWeights :: Strategy -> PayoffVector -> VU.Vector Double
baseWeights RandomLearning payoffVector = VU.replicate (VU.length payoffVector) 1.0
baseWeights PayoffBasedLearning payoffVector = payoffVector

-- | Computes the weight vector w(s, r) based on strategy and current repertoire
weightVector :: Strategy -> Repertoire -> PayoffVector -> VU.Vector Double
weightVector strategy repertoire payoffVector =
    let w_star = baseWeights strategy payoffVector
        unlearned = VU.map not repertoire
        w_unlearned = VU.zipWith (\w u -> if u then w else 0.0) w_star unlearned
        total = VU.sum w_unlearned
    in if total == 0
        then VU.replicate (VU.length repertoire) 0.0
        else VU.map (/ total) w_unlearned

-- | Computes transition probabilities from the given repertoire
transitionFromState :: Strategy -> Repertoire -> AdjacencyMatrix -> PayoffVector -> [(Repertoire, Double)]
transitionFromState strategy repertoire adjMatrix payoffVector =
    let w = weightVector strategy repertoire payoffVector
        learnable = learnability repertoire adjMatrix
        transitions =
            [ (learnTrait repertoire j, w VU.! j)
            | j <- [0 .. VU.length repertoire -1]
            , learnable VU.! j
            , let wj = w VU.! j, wj > 0.0
            ]
    in transitions
  where
    learnTrait :: Repertoire -> Int -> Repertoire
    learnTrait r j = r VU.// [(j, True)]

-- | Computes probability of staying in the same repertoire
stayProbability :: Strategy -> Repertoire -> AdjacencyMatrix -> PayoffVector -> Double
stayProbability strategy repertoire adjMatrix payoffVector =
    let transitions = transitionFromState strategy repertoire adjMatrix payoffVector
        totalTransitionProb = sum (map snd transitions)
    in 1.0 - totalTransitionProb

-- | Generates all reachable repertoires starting from the initial repertoire
generateReachableRepertoires :: Strategy -> AdjacencyMatrix -> PayoffVector -> [Repertoire]
generateReachableRepertoires strategy adjMatrix payoffs = go [initialRepertoire] []
  where
    n = V.length adjMatrix
    initialRepertoire = VU.accum (const id) (VU.replicate n False) [(0, True)]  -- Root trait is learned
    go [] visited = visited
    go (r:queue) visited
      | r `elem` visited = go queue visited  -- Skip if already visited
      | otherwise =
          let transitions = [ r' | (r', _) <- transitionFromState strategy r adjMatrix payoffs ]
              newQueue = queue ++ transitions
          in go newQueue (visited ++ [r])
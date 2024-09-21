{-# LANGUAGE FlexibleContexts #-}
module Payoff (
    generatePayoffs
) where

import qualified Data.Vector.Unboxed as VU
import Control.Monad.State.Strict (runState, state, State)
import System.Random (StdGen, randomR)
import Types (PayoffVector)

-- | Generates payoffs based on distances and alpha parameter using a random seed
generatePayoffs :: VU.Vector Int -> Double -> StdGen -> (PayoffVector, StdGen)
generatePayoffs distances alpha gen =
    let (payoffsList, gen') = runState (mapM generatePayoff' $ VU.toList distances) gen
    in (VU.fromList payoffsList, gen')
  where
    generatePayoff' :: Int -> State StdGen Double
    generatePayoff' dist
        | dist == 0 = return 0.0 -- Payoff of the root trait is always 0
        | otherwise = do
            randValue <- state $ randomR (0.0, 1.0)
            let payoff = randValue + alpha * fromIntegral dist
            return payoff
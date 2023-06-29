
# Problem 2 ----
## 2.1

strategy1 <- function(ownpastactions, partnerpastactions) { return(TRUE) } # always True

strategy2 <- function(ownpastactions, partnerpastactions) { return(runif(1) > 0.75) } # 3/4 False

strategy3 <- function(ownpastactions, partnerpastactions) { # 1/3 True
    n <- length(ownpastactions)
    
    if (n %% 3 == 0) {
        return(TRUE) 
    } else { 
        return(FALSE) 
    } 
}

## 2.2

payout <- function(action_p1, action_p2, payoutvector = c(5, 3, 1, 0)) {
    # best (5) > good (3) > bad (1) > worst (0)
    if (action_p1 & action_p2) return(c(payoutvector[2], payoutvector[2])) 
    if (action_p1 & ! action_p2) return(c(payoutvector[4], payoutvector[1])) 
    if (! action_p1 & action_p2) return(c(payoutvector[1], payoutvector[4])) 
    if (! action_p1 & ! action_p2) return(c(payoutvector[3], payoutvector[3]))
}

simulate_2P_PD <- function(strategy_p1, strategy_p2, N, payoutvector = c(5, 3, 1, 0)) { 
    total_payout <- c(0, 0)
    past_actions_p1 <- vector()
    past_actions_p2 <- vector()
    for (i in 1:N) { 
      
      out1 <- strategy_p1(past_actions_p1, past_actions_p2)
      out2 <- strategy_p2(past_actions_p2, past_actions_p1)
      
      local_payout <- payout(out1, out2)
      
      total_payout <- total_payout + local_payout
      
      past_actions_p1 <- c(past_actions_p1, out1)
      past_actions_p2 <- c(past_actions_p2, out2)
        
    }
    return(total_payout) 
}

## 2.3

tournament <- function(list_of_strategies = list(strategy1, strategy2, strategy3), 
                       N = 50, Nsim = 100, payoutvector = c(5, 3, 1, 0), 
                       PLAY_YOURSELF = FALSE) {
    
    n <- length(list_of_strategies) 
    total_payouts = rep(0, n) 
    av_payout_mat = matrix(0, n, n)
    
    for (i in 1:(n - 1 + PLAY_YOURSELF)) {
        for (j in (i + 1 - PLAY_YOURSELF):n) {
            average_payout <- ?? #Remember to average 100 repetitions 
                av_payout_mat[i,j] <- ??
                    av_payout_mat[j,i] <- ??
                        total_payouts[i] <- total_payouts[i] + average_payout[1] / 
                            (n - 1 + PLAY_YOURSELF) 
                        total_payouts[j] <- total_payouts[j] + average_payout[2] / 
                            (n - 1 + PLAY_YOURSELF)
        } 
    }
    return(av_payout_mat) 
}

heatmap(av_payout_mat, scale = 'none', revC = TRUE, Rowv = NA, Colv = NA, margins = c(10, 10))


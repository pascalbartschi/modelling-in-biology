import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression


def payout(action_p1, action_p2, pay = np.array([5, 3, 1, 0])):

    if action_p1 and action_p2:
        return np.array([pay[1], pay[1]])

    if action_p1 and not action_p2:
        return np.array([pay[3], pay[0]])

    if not action_p1 and action_p2:
        return np.array([pay[0], pay[3]])

    if not action_p1 and not action_p2:
        return np.array([pay[2], pay[2]])



def simulate_2P_PD(strategy_p1, strategy_p2, N, pay = np.array([5, 3, 1, 0])):

    total_pay = np.array([0, 0])

    past_actions_p1 = np.array([])
    past_actions_p2 = np.array([])

    for i in range(N):

        out1 = strategy_p1(past_actions_p1, past_actions_p2)
        out2 = strategy_p2(past_actions_p2, past_actions_p1)

        total_pay += payout(out1, out2)

        past_actions_p1 = np.append(past_actions_p1, out1)
        past_actions_p2 = np.append(past_actions_p2, out2)

    return total_pay

def tournament(list_of_strategies, N, Nsim, pay = np.array([5, 3, 1, 0])):

    n = len(list_of_strategies)

    total_payouts = np.zeros(n)
    av_payout_mat = np.zeros((n ,n))

    for i in range(n):
        for j in range(i+1, n):
            sum_payout = np.array([0, 0])
            for sim in range(Nsim):
                sum_payout += simulate_2P_PD(list_of_strategies[i], list_of_strategies[j], N)

            # avg performance per battle
            average_payout = sum_payout / Nsim
            total_payouts[i] += average_payout[0]
            total_payouts[j] += average_payout[1]

            # append to matrix
            av_payout_mat[i, j] = average_payout[0]
            av_payout_mat[j, i] = average_payout[1]


    return total_payouts, av_payout_mat



def tit_for_tat(ownpast, partnerpast):

    n = len(ownpast)

    if n == 0:
        return True

    else:
        return partnerpast[-1]


def tit_for_two_tats(ownpast, partnerpast):

    n = len(ownpast)

    if n == 0:
        return True

    if not partnerpast[-2:].any():
        return False

    else:
        return True

def two_tits_for_tat(ownpast, partnerpast):

    n = len(ownpast)

    if n == 0:
        return True

    if partnerpast[-2:].all():
        return False

    else:
        return True











import sys
from loguru import logger
import numpy as np
import pandas as pd
import random
from scipy.stats import norm

# logger.remove(0)

logger.add(
    'apamix.log',
    rotation='10 MB',
    colorize=True,
    level="DEBUG"
)

logger.add(sys.stderr, level="DEBUG", colorize=True)

"""
This was the prototype of apa mixture module

required args for running this module.

Known variable:
    - Start site of R2, denote as x_n
    - Mapped length of R2, denote as l_n
    - Length of R1 which was poly(A), denote as r_n

Latent variable:
    - The real poly(A) length at the capture site, denote as s_n
    - The given reads were from the specific poly(A), denote as sznk

"""


def lik_l_x(l, x, theta):
    """

    For given l_n and x_n, return the likelihood of theta_k

    """
    utr_len = theta - x + 1  # given poly(A) site substract the start site infer to the UTR length

    valid_inds = l <= utr_len  # The reads must at the 5'

    res = (valid_inds).astype('float64')

    if any(np.isnan(res)):
        logger.warning("Some reads' length were 0, please check it.\n Location at calculating the L(l_n|x_n, theta_k)")

    res[valid_inds] = 1 / utr_len[valid_inds]

    if any(np.isinf(res)):
        logger.warning("Results contain Inf value.\n Location at calculating the L(l_n|x_n, theta_k)")
        sys.exit(1)

    return res


def lik_x_s(x, s, theta, mu_f, sigma_f):
    """
    For given x_n and s_n
    """

    res = norm.pdf(x,
                   theta + s - mu_f + 1,
                   sigma_f
                   )

    if any(np.isinf(res)):
        logger.warning("Results contain Inf value.\n Location at calculating the L(x_n|s_n, theta_k)")
        sys.exit(1)

    return res


def lik_r_s(r, s):
    """
    poly(A) length
    """

    res = (r <= s) / s

    if any(np.isinf(res)):
        logger.warning("Results contain Inf value.\n Location at calculating the L(r_n|s_n)")
        sys.exit(1)

    return res


def lik_s(s_dis):
    """
    uniform distribution for the poly(A) length
    """
    n = len(s_dis)

    res = np.repeat(1 / n, n)

    if any(np.isinf(res)):
        logger.warning("Results contain Inf value.\n Location at calculating the L(r_n|s_n)")
        sys.exit(1)

    return res


def likelihood(x_arr, l_arr, r_arr, s_arr, theta, s_dis_arr, mu_f, sigma_f):
    """
    s_dis_arr = np.arange(start = 20,stop = 250 + 1, step = 10)

    get the all likelihood results, and calculate the expectation.
    1. p(l_n|x_n, theta_k) -> lik_l_x
    2. p(x_n|s_n,theta_k) -> lik_x_s
    3. p(r_n|s_) -> lik_r_s
    """

    res = np.repeat(0.0, len(x_arr))

    ps = lik_s(s_dis_arr)
    ns = len(s_dis_arr)

    # paired alignment reads
    valid_inds = ~np.isnan(s_arr)
    oth_inds = ~valid_inds

    lik_l_x_ = lik_l_x(l=l_arr[valid_inds],
                       x=x_arr[valid_inds],
                       theta=theta)

    lik_x_s_ = lik_x_s(x=x_arr[valid_inds],
                       s=s_arr[valid_inds],
                       theta=theta,
                       mu_f=mu_f,
                       sigma_f=sigma_f)

    res[valid_inds] = lik_l_x_ * lik_x_s_

    if not any(oth_inds):
        return res
    for i in range(ns):
        lik_l_x_ = lik_l_x(l=l_arr[oth_inds],
                           x=x_arr[oth_inds],
                           theta=theta)

        lik_x_s_ = lik_x_s(x=x_arr[oth_inds],
                           s=s_dis_arr[i],
                           theta=theta,
                           mu_f=mu_f,
                           sigma_f=sigma_f)

        lik_r_s_ = lik_r_s(r=r_arr[oth_inds],
                           s=s_dis_arr[i])

        res[oth_inds] = res[oth_inds] + lik_l_x_ * lik_x_s_ * lik_r_s_ * ps[i]

        if np.isnan(sum(res)):
            logger.error("The res contains NaN.")
            sys.exit(1)
        if any(np.isinf(res)):
            logger.error("The res contains Inf.")
            sys.exit(1)

    return (res)


def cal_theta(Z, ws, x_arr, l_arr, r_arr, s_arr, theta_arr, k, L):
    """
    k max more than zero
    """

    k_ = len(ws)

    min_theta = min(theta_arr) if k == 0 else theta_arr[k - 1] + 1
    max_theta = L if k == k_ - 1 else theta_arr[k] - 1

    all_theta_arr = np.arange(min_theta, max_theta, 10)

    pa_gap = 9

    all_theta_arr = all_theta_arr[
        np.logical_and(np.abs(all_theta_arr - min_theta) > 9,
                       np.abs(all_theta_arr - max_theta) > 9)]

    pass


def e_step(ws, x_arr, l_arr, r_arr, s_arr, theta_arr, s_dis_arr, sigma_f, mu_f):
    """

    Here need to reformat the name of variable

    """

    k_ = len(ws)  # numbers of poly(A) site
    n_ = len(x_arr)  # numbers of reads
    z_ = np.zeros((n_, k_))

    for k in range(k_):
        z_[:, k] = ws[k] * likelihood(x_arr=x_arr,
                                      l_arr=l_arr,
                                      r_arr=r_arr,
                                      s_arr=s_arr,
                                      theta=theta_arr[k],
                                      s_dis_arr=s_dis_arr,
                                      mu_f=mu_f,
                                      sigma_f=sigma_f)

    res = z_ / z_.sum(1)[:, None]

    return (res)


def m_step(Z, ws, x_arr, l_arr, r_arr, s_arr, theta_arr, k):
    """
    Calulation the ws

    """
    k_ = len(x_arr)
    n_ = len(ws)

    new_ws = Z.sum(0) / k_
    new_thera_arr = theta_arr.copy()
    if k > 0:
        new_thera_arr[k] = cal_theta(Z, ws, x_arr, l_arr, r_arr, s_arr, theta_arr, k)

    return [new_ws, new_thera_arr]


def elbo(Z, ws, x_arr, l_arr, r_arr, s_arr, theta_arr):
    """
    calculating entropy
    """

    Z_copy = Z.copy()
    Z_copy[Z != 0] = np.log(Z[Z != 0])

    entropy = -1 * Z * Z_copy

    elb = exp_log_lik(Z, ws, x_arr, l_arr, r_arr, s_arr, theta_arr) + entropy.sum()

    if np.isnan(elb):
        logger.error("lower bounder is Na.")
        sys.exit(0)

    return (elb)


def exp_log_lik(Z, ws, x_arr, l_arr, r_arr, s_arr, theta_arr):
    if np.isnan(theta_arr).any():
        logger.warning("Some thera is NaN.")

    k_ = len(ws)  # numbers of poly(A) site
    n_ = len(x_arr)  # numbers of reads
    z_ = np.zeros((n_, k_))
    for k in range(k_):
        z_[:, k] = np.log(ws[k] * likelihood(x_arr=x_arr,
                                             l_arr=l_arr,
                                             r_arr=r_arr,
                                             s_arr=s_arr,
                                             theta=theta_arr[k],
                                             s_dis_arr=s_dis_arr,
                                             mu_f=mu_f,
                                             sigma_f=sigma_f))
    z_[Z == 0] = 0
    ell = z_ * Z
    if np.isnan(ell).any():
        logger.error("Exp log likelihood contain Nan.")
        sys.exit(1)

    return (ell.sum())


def gen_k_arr(K, n):
    """
    generate non same elements
    """

    def random_sel(K, trial=200):
        count_index = 0
        pool = np.arange(K)
        last = None
        while count_index < trial:
            count_index += 1
            random.shuffle(pool)
            if pool[0] == last:
                swap_with = random.randrange(1, len(pool))
                pool[0], pool[swap_with] = pool[swap_with], pool[0]
            for item in pool:
                yield item

            last = pool[-1]

    if K <= 1:
        return np.repeat(K, n)
    else:
        k_lst = list(random_sel(K, trial=n))
        return np.array(k_lst)


def init_ws(n_apa):
    """

    number of apa site

    """

    ws = np.random.rand(n_apa) + 1
    ws = ws / ws.sum()
    return (ws)


def init_thera(x_arr, l_arr, L, K, pre_theta_arr=None, n_max_trial=200, min_gap_pa_sites=50):
    """
    """
    min_theta = min(l_arr)
    min_max_theta = max(x_arr + l_arr + 1)
    max_theta = L
    full_theta = np.arange(start=min_theta, stop=max_theta + 1, step=10)
    if pre_theta_arr:
        logger.error("For now, not support the pre_theta_arr mode!")
        sys.exit(1)
    else:
        tmp_k = 0
    for trial in range(n_max_trial):
        if not pre_theta_arr:
            tmp_theta_arr = np.random.choice(full_theta, K - tmp_k)
        else:
            tmp_theta_arr = np.append(np.random.choice(full_theta, K - tmp_k), pre_theta_arr)
        tmp_theta_arr.sort()
        flag1 = all(np.diff(tmp_theta_arr) >= min_gap_pa_sites)
        flag2 = tmp_theta_arr[K - 1] >= min_max_theta
        if flag1 & flag2:
            break
    if trial == n_max_trial + 1:
        logger.error("Failed to generate valid thera afer {} trials".format(n_max_trial + 1))
        sys.exit(1)
    else:
        pass
    return tmp_theta_arr


def main():
    file = 'a'
    df = pd.read_table(file)
    n_max_apa = 3
    n_min_apa = 1

    x_arr = df.V1.values
    l_arr = df.V2.values
    r_arr = df.V3.values
    s_arr = df.V6.values

    L = 1100

    min_ws = 0.1
    min_gap_pa_sites = 50
    pre_theta_arr = None
    n_apa = 3

    min_la = 20
    max_la = 250

    s_dis_arr = np.arange(start=min_la, stop=max_la + 1, step=10)
    k_arr = gen_k_arr(len(theta_arr), nround)
    nround = 20

    mu_f = 300
    sigma_f = 50

    ## variable from init
    ws = init_ws(n_apa)
    theta_arr = init_thera(x_arr, l_arr, L, n_apa)

    Z = e_step(ws=ws,
               x_arr=x_arr,
               l_arr=l_arr,
               r_arr=r_arr,
               s_arr=s_arr,
               theta_arr=theta_arr,
               s_dis_arr=s_dis_arr,
               mu_f=mu_f,
               sigma_f=sigma_f)

    lb_new = elbo(Z, ws, x_arr, l_arr, r_arr, s_arr, theta_arr)
    # res = m_step(Z, ws, x_arr, l_arr, r_arr, s_arr,  theta_arr, k_arr[i])
    res = m_step(Z, ws, x_arr, l_arr, r_arr, s_arr, theta_arr, 0)


if __name__ == '__main__':
    main()




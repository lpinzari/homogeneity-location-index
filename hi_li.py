import math
from typing import List

# Title       : Python implementation of Ludovico Pinzaris Homogeneity index (HI) and Location Index (LI)
# Author      : Linus Kohl <linus@munichresearch.com>
# Python      : > 3.6


def concI(p: List[int]) -> float:
    """
    Concentration Index (CI)
    Function to compute the concentration index (CI) for a vector of observations
    (e.g. the number of people in each quantile of the socioeconomic index).
    Args:
        p(List[int]): Numeric vector with non-negative integer values
    Returns:
        float: The Concentration Index (CI) of p. (real number [0 1])
    Example:
        uniform distribution
        p = [10,10,10,10,10,10,10,10,10,10]
        r = concI(p) # 0
        p = [50,0,0,0,0,0,0,0,0,0,50]
        r = concI(p) # 0.88
        p = [100,0,0,0,0,0,0,0,0,0]
        r = concI(p) # 1
    """
    d = len(p)
    pop = sum(p)
    p = [x / pop for x in p]
    p.sort()  # sort the pdf vector in ascending order
    # Compute the cumulative frequencies Lorenz Curve
    lc = [0] * (d + 1)  # init array
    lc[1] = lc[0] + p[0]
    for i in range(1, d):
        lc[i + 1] = lc[i] + p[i]
    # Compute the Area vector Under the Lorenz Curve
    b = 1 / d  # base of the trapezoid
    imp_area = 1 - b  # to normalize the result
    la = [0] * d  # init array
    for i in range(0, d):
        la[i] = b * (lc[i] + lc[i + 1]) / (2 * imp_area)
    # Compute the Lorenz Curve of the unifotm distr.
    pu = [1 / d] * d  # init
    lcu = [0] * (d + 1)  # init
    lcu[1] = lcu[0] + pu[0]
    for i in range(1, d):
        lcu[i + 1] = lcu[i] + pu[i]
    # Compute the area under the Lorenz Curve of the unifoirm distribution
    lau = [0] * d
    for i in range(0, d):
        lau[i] = b * (lcu[i] + lcu[i + 1]) / (2 * imp_area)
    # Compute the Zonoid Area (Area between the lorenz Curve and Uniform distribution)
    area = 0
    for i in range(0, d):
        area = area + lau[i] - la[i]
    # rounding the value
    area = round(area, 4)

    return 2 * area


def conv(x: List[float], y: List[float]) -> List[float]:
    """
    Unidimensional Convolution
    Function to compute the convolution of two vectors
    Args:
        x(List[int]): numeric vector representing polynomial coefficient (1 x m)
        y(List[int]): numeric vector representing polynomial coefficient (1 x n)
    Returns:
       List[float]: Coefficient vector resulting from multiplying the polynomial represented
                    by x by the polynomial represented by y (1 x m+n-1)
    Example:
          x = [0.25,0.25,0.25,0.25]
          y = [1,1,1,1]
          r = conv(x,y) # 0.25 0.50 0.75 1.00 0.75 0.50 0.25
          x = [1,0,0,0]
          y = [1,1,1,1]
          r = conv(x,y) # 1 1 1 1 0 0 0
    """
    m = len(x)
    n = len(y)
    z = [0] * (m + n - 1)
    for j in range(0, m):
        for k in range(0, n):
            z[j + k] = z[j + k] + x[j] * y[k]

    return (z)


def corr(x: List[float]) -> List[float]:
    """
    Unidimensional autocorrelation
    Uses conv to compute the autocorrelation of a vector
    Params:
        x(List[int]): Numeric vector representing polynomial coefficient (1 x n)
    Returns:
        List[float]: The coefficient vector resulting from the autocorrelation (1 x 2n-1)
    Example:
        x = [0.25,0.25,0.25,0.25]
        r = corr(x) # 0.0625 0.1250 0.1875 0.2500 0.1875 0.1250 0.0625
    """
    R = conv(x, x[::-1])

    return R


def div(x: List[float]) -> float:
    """
    Divergence Index
    Uses conv and corr to compute the polarization divergence of a
    probability vector.
    Params:
        x(List[float]): Numeric vector representing polynomial coefficient
                        distribution (pdf)
    Returns:
        float: The Divergence index. (1 x 1) real number [0 1]
    Examples:
        x = [0.5,0,0,0,0,0,0,0,0,0.5]
        r = div(x) # 0.2973122
        x = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1] # Uniform distribution
        r = div(x) # 0.1229451
        x = [0,0,0,1,0,0,0,0,0,0]
        r = div(x) # 0
    """
    d = len(x)
    comb = [1] * d
    # compute the Bilateral Cumulative Distributive function and its autocorrelation (spectra)
    bcdf = conv(x, comb)
    sbcdf = corr(bcdf)
    # compute the Singleton autocorrelation function
    imp = [0] * d
    imp[0] = 1
    bcdfi = conv(imp, comb)
    sbcdfi = corr(bcdfi)
    # normalized energy (pdf of the signal)
    E = sum(sbcdfi)
    sbcdf = [x / E for x in sbcdf]
    sbcdfi = [x / E for x in sbcdfi]
    # compute the binary logarithm of the signals
    l = len(sbcdf)
    lgs1 = [0] * l
    lgsi = [0] * l
    for i in range(0, l):
        if sbcdf[i] > 0:
            lgs1[i] = math.log(sbcdf[i], 2)
        if sbcdfi[i] > 0:
            lgsi[i] = math.log(sbcdfi[i], 2)
    # compute the -log M(X)
    M = [0] * l
    for i in range(0, l):
        M[i] = (sbcdf[i] + sbcdfi[i]) / 2
        if M[i] > 0:
            M[i] = -math.log(M[i], 2)
    # compute D(I,M) D(I,S)
    div = 0
    for i in range(0, l):
        i_m = sbcdfi[i] * (lgsi[i] + M[i])
        i_s = sbcdf[i] * (lgs1[i] + M[i])
        div = div + i_m + i_s

    return div


def divConst(n: int) -> float:
    """
    Divergence Index Constant
    Uses div to compute the distribution divergence constant for the maximum
    variance.
    Params:
        n(int): The number of bins in the pdf (pdf)
    Returns:
        float: The Divergence index constant. (1 x 1) real number [0 1)
    Examples:
        x = [0.5,0,0,0,0,0,0,0,0,0.5]
        r = div(x)       # 0.2973122
        r = divConst(10) # 0.2973122
    """
    # create a bimodal distribution
    p = [0] * n  # pdf
    p[0] = 0.5
    p[n - 1] = 0.5
    # compute the divergence index of the bimodal pdf
    const = div(p)

    return const


def hom(p: List[int]) -> float:
    """
    Homogenity Index
    Uses concI and div to compute the Homogeneity Index for a vector of
    observations.
    Params:
        p(List[int]): Numeric vector with non-negative integer values. (1 x n)
    Returns:
        float: The Homogeneity Index (HI) of p: real number  [0 1]
    Examples:
        p = [10,10,10,10,10]
        r = hom(p) # 0
        p = [0,0,50,0,0]
        r = hom(p) # 1
    """
    # compute the concentration index
    conc = concI(p)
    # compute the pdf of p
    d = len(p)
    pop = sum(p)
    p = [x / pop for x in p]
    # compute the Divergence Index of the distribution
    div_i = div(p)
    # compute the Divergence Index for the Uniform distribution
    uniform = [1 / d] * d
    divUnif = div(uniform)
    # compute the Homogeneity Index
    H = (conc + divUnif - div_i) / (1 + divUnif)

    return H


def hom_mn(m: int, n: int) -> float:
    """
    Homogeneity Index - True diversity
    Uses hom to compute the Homogenity Index for a distribution of m bins and
    n equally abundant contiguous categories: pdf_mn
    Params:
        m(Int): Integer indicating the number of bins in the distribution.
        n(Int): integer indicating the number of contiguous bins with equal
                number of observations
    Returns:
        float: The Homogeneity Index (HI) of pdf_mn: real number  [0 1]
    Examples:
        r = hom_mn(5,5) # 0
        r = hom_mn(5,1) # 1
        p = [0.5,0.5,0,0,0]
        r = hom(p)   # 0.76
        r = hom(5,2) # 0.76
    """
    ## create a pdf_mn
    p = [0] * m
    for j in range(0, n):
        p[j] = 1 / n
    print(p)
    h = hom(p)

    return h


def hom_class(a: float, b: float, c: float, p) -> str:
    """
    Homogeneity Index - Classification
    Uses hom to compute the Homogenity class for a distribution
    Params:
        a(float): Real number indicating the lower bound of the first class - A
        b(float): Real number indicating the lower bound of the second class - B
        c(float): Real number indicating the lower bound of the third class - C
    Returns:
        str: Homogenity Class
    """
    h = hom(p)
    if h >= a:
        rclass = 'A'
    elif h >= b:
        rclass = 'B'
    elif h >= c:
        rclass = 'C'
    else:
        rclass = 'D'

    return (rclass)


def locVec(x: List[float]) -> List[float]:
    """
    Location Index Supporting Function
    Function to compute the Concentration Location Index vector of a pdf
    Params:
        x(List[float]):  Probability density function vector. (1 x n)
    Returns:
        List[float]: The Location Index vector score (LIS) of x
    Examples:
        x = [0.5,0,0,0,0.5]
        v = locVec(x) # 0.6 0.6 0.6 0.6 0.6
        x = [1,0,0,0,0]
        x = locVec(x) # 1 0.8 0.6 0.4 0.2
    """
    n = len(x)
    k = [0] * n  # bins scores
    for i in range(0, n):
        j = 0  # iterator for the nested intervals (width)
        s = 0
        vs = [0] * n
        for j in range(0, n):
            if j == 0:
                s = x[i]  # interval width zero (initial point)
            else:
                fw = i + j  # interval border right
                bk = i - j  # interval border left
                if bk >= 0 and fw < n:
                    s = s + x[bk] + x[fw]
                elif bk >= 0:
                    s = s + x[bk]
                elif fw < n:
                    s = s + x[fw]
            vs[j] = s
        k[i] = sum(vs)
    z = 1 / n
    k = [x * z for x in k]

    return k


def loc(x):
    """
    Location Index
    Function to compute the Location Index (LI) and Compactness(C).
    LI gives the minimum and maximum position of the bins with the maximum concentration
    Params:
        x(List[float]): Probability density function vector. (1 x n)
    Returns:
        List[int]: The Location Index (LI) and Compactness of x (1x2)
    Examples:
        x = c(0.5,0,0,0,0.5]
        v = loc(x) # 1 5
        x = c(1,0,0,0,0]
        x = loc(x) # 1 1
    """
    n = len(x)
    loc1 = 0  # minimum position of the bin with maximum score
    loc2 = 0  # maximum position of the bin with maximum score
    ## compute the Location Index vector and the maximum score
    v = locVec(x)
    m = max(v)
    for j in range(0, n):
        if v[j] == m:
            if loc1 > 0:
                loc2 = j
            else:
                loc1 = j
                loc2 = j

    cl = [loc1, loc2]

    return cl

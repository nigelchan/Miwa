Miwa
====

Functional Programming in Finance

Demonstration:
First calibrate a heston model to market data, then use it to price a European call option.

    [DATA] S: 4468.17
    [DATA] Implied Volatility Surface:
    K\ T         0.0356         0.1123         0.2055         0.4521         0.7014         0.9452         1.4356         1.9260
    ============================================================================================================================
    3400         0.6625         0.4875         0.4204         0.3667         0.3431         0.3267         0.3121         0.3121
    3600         0.6007         0.4543         0.3967         0.3511         0.3279         0.3154         0.2984         0.2921
    3800         0.5084         0.4221         0.3718         0.3327         0.3155         0.3027         0.2919         0.2889
    4000         0.4541         0.3869         0.3492         0.3149         0.2963         0.2926         0.2819         0.2800
    4200         0.4060         0.3607         0.3330         0.2999         0.2887         0.2811         0.2751         0.2775
    4400         0.3726         0.3396         0.3108         0.2781         0.2788         0.2722         0.2661         0.2686
    4500         0.3550         0.3277         0.3012         0.2781         0.2781         0.2661         0.2661         0.2681
    4600         0.3428         0.3209         0.2958         0.2740         0.2688         0.2627         0.2580         0.2620
    4800         0.3302         0.3062         0.2799         0.2631         0.2573         0.2533         0.2504         0.2544
    5000         0.3343         0.2959         0.2705         0.2540         0.2504         0.2464         0.2448         0.2462
    5200         0.3460         0.2845         0.2624         0.2463         0.2425         0.2385         0.2373         0.2422
    5400         0.3857         0.2860         0.2578         0.2399         0.2357         0.2327         0.2312         0.2351
    5600         0.3976         0.2860         0.2607         0.2356         0.2297         0.2268         0.2241         0.2320

    Start calibration...
                            Theta           Kappa           Sigma             Rho              v0
    iter:   1 x =      0.06902680      1.95039631      1.11171680     -0.48311396      0.12038514
    iter:   2 x =      0.07965551      5.89670250      2.00306657     -0.51826382      0.15376352
    iter:   3 x =      0.07316429     10.60620585      2.70947333     -0.51846056      0.17536640
    iter:   4 x =      0.07459079     13.60843403      3.07934324     -0.51511396      0.18614506
    iter:   5 x =      0.07474602     14.96884193      3.23297170     -0.51336425      0.18961481
    iter:   6 x =      0.07481433     15.43028832      3.28430133     -0.51265976      0.19067708
    iter:   7 x =      0.07483687     15.56804037      3.29963896     -0.51243064      0.19098784
    iter:   8 x =      0.07484359     15.60736268      3.30402290     -0.51236303      0.19107631
    iter:   9 x =      0.07484551     15.61844126      3.30525838     -0.51234380      0.19110123
    iter:  10 x =      0.07484605     15.62154977      3.30560507     -0.51233838      0.19110822
    iter:  11 x =      0.07484620     15.62241905      3.30570212     -0.51233686      0.19111018

    status = converged


    Option pricing (K = 4200; T = 0.4583)...

    Using Crude Monte Carlo Engine...
    Call price = 533.834844

    Using Heston Semi Analytic Engine...
    Call price = 531.776813

    Using Heston Monte Carlo Engine...
    Call price = 528.440436

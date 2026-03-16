
Old Source: randomPointTrapOptimumCleanable

Old Tracker: fixedEffDaggerHitTime

Output that matters:
1. T_start
2. T_final
3. E_start
4. E_final
5. X_init
6. Y_init
7. Z_init
8. Vx_init
9. Vy_init
10. Vz_init
11. X_final
12. Y_final
13. Z_final
14. Vx_final
15. Vy_final
16. Vz_final
17. zOff
18. nHit
19. nHitHouseLow
20. nHitHouseHigh
21. deathTime

v1 baseline:
 - source: Randomized Source
 - field: HalbachField
 - integrator: SymplecticIntegrator
 - tracker: FixedEffDaggerTracker
 - writer: CsvWriter
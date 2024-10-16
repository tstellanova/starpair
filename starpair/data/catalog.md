## From `download_oneshot.py`:

### Limit 2e8 r45 d100
-  antigalactic line stars written to:
`antigalactic_L2e+08_r45_d100_18000_0.fits.gz`
- 26742 nearby galactic line stars written to:
`galactic_L2e+08_r45_d100_0_0.fits.gz`
query >>> elapsed: 963.44 seconds

### Limit 2e9 r30 d100 
- 23043 nearby antigalactic line  stars 23043 written to:
`antigalactic_L2e+09_r30_d100_18000_0.fits.gz`
- 185819 nearby galactic line stars 185819 written to:
`galactic_L2e+09_r30_d100_0_0.fits.gz`

### Limit 2e9 r45 d100
- 48954 nearby antigalactic line stars written to:
`antigalactic_L2e+09_r45_d100_18000_0.fits.gz`
antigalactic query >>> elapsed: 2031.90 seconds
- 240634 nearby galactic line stars written to:
`galactic_L2e+09_r45_d100_0_0.fits.gz`
galactic query >>> elapsed: 3689.79 seconds

### Limit 2e9 outer cone ro60 ri45
- 72071 nearby stars  written to:
`galactic_outercone_L2e+09_ri45_ro60_d100_0_0.fits.gz`
galactic outercone query >>> elapsed: 6893.27 seconds
### Limit 2e9 outer cone ro90 ri60
- 97095 nearby stars  written to:
`galactic_outercone_L2e+09_ri60_ro90_d100_0_0.fits.gz`
  galactic outercone query >>> elapsed: 2078.54 seconds

## Re-query oneshot with L,B parameters 45 degrees
DESIGNATION,ref_epoch,ra,dec,parallax,l,b,pm,pmra,pmdec,phot_g_mean_mag, ABS(1000./parallax) AS dist_pc, DISTANCE(180.00, 0.00, l, b) AS ang_sep

### Limit 2e9 igalactic r45
240634 nearby stars  written to:
`galactic_L2e+09_r45_d100_0_0.fits.gz`
galactic query >>> elapsed: 2274.69 seconds

### Limit 2e9 antigalactic r45 
48954 nearby stars  written to:
`antigalactic_L2e+09_r45_d100_18000_0.fits.gz`
antigalactic query >>> elapsed: 1817.21 seconds


## Re-query oneshot with L, B parameters 

DESIGNATION,ref_epoch,ra,dec,parallax,l,b,pm,pmra,pmdec,phot_g_mean_mag, ABS(1000./parallax) AS dist_pc, DISTANCE(180.00, 0.00, l, b) AS ang_sep

### Limit 2e9 antigalactic ri45 ro60
34929 nearby stars  written to:
`antigalactic_outercone_L2e+09_ri45_ro60_d100_18000_0.fits.gz`
galactic outercone query >>> elapsed: 3756.50 seconds

### Limit 2e9 antigalactic ri60 ro90
80848 nearby stars  written to:
`antigalactic_outercone_L2e+09_ri60_ro90_d100_18000_0.fits.gz`
galactic outercone query >>> elapsed: 3236.75 seconds

### Limit 2e9 galactic ri45 ro60
72071 nearby stars  written to:
`galactic_outercone_L2e+09_ri45_ro60_d100_0_0.fits.gz`
galactic outercone query >>> elapsed: 1829.84 seconds

### Limit 2e9 galactic ri60 ro90
97095 nearby stars  written to:
`galactic_outercone_L2e+09_ri60_ro90_d100_0_0.fits.gz`
galactic outercone query >>> elapsed: 1847.82 seconds

### Limit 2e9 galactic r90 d100 
409800 nearby stars written to:
`galactic_L2e+09_r90_d100_0_0.fits.gz`
galactic query >>> elapsed: 3595.14 seconds

### Limit 2e9 antigalactic r90 d100
164731 nearby stars written to:
`antigalactic_L2e+09_r90_d100_18000_0.fits.gz`
antigalactic query >>> elapsed: 6189.96 seconds

## Approximate quantities of data for 30, 45, 90 degree cones
```shell
(astrodata_env) todd@relap data % ls -latr galactic_L2e+09_*d100_0_0.fits.gz
-rw-r--r--  1 todd  staff  13667577 Oct  3 20:08 galactic_L2e+09_r30_d100_0_0.fits.gz
-rw-r--r--  1 todd  staff  21606955 Oct  4 23:21 galactic_L2e+09_r45_d100_0_0.fits.gz
-rw-r--r--  1 todd  staff  37016369 Oct  5 01:42 galactic_L2e+09_r90_d100_0_0.fits.gz
```
```shell
(astrodata_env) todd@relap data % ls -latr antigalactic_L2e+09_r*d100_18000_0.fits.gz
-rw-r--r--  1 todd  staff   1734133 Oct  3 19:38 antigalactic_L2e+09_r30_d100_18000_0.fits.gz
-rw-r--r--  1 todd  staff   4458179 Oct  4 23:51 antigalactic_L2e+09_r45_d100_18000_0.fits.gz
-rw-r--r--  1 todd  staff  15010706 Oct  5 10:28 antigalactic_L2e+09_r90_d100_18000_0.fits.gz
```


Assume my CSV file has the following comma-separated fields on each line:
max_node_dist: float64
min_node_dist: float64
ang_sep: float64
radial_sep: float64
src_coord_str: two space-separated float64 numbers in a string of the form '175.521 -75.6869' representing g_lon and g_lat float64 fields
dst_coord_str: same format as src_coord_str
src_id: uint64
dst_id: uint64




# Example: Repeat masking chromosome 28 of the Naked Mole Rat with Earl Grey

## Installing Earl Grey

Create a conda environment for Earl Grey

```
conda create -n earlgrey_env earlgrey=4.4.0
```

Activate conda environment

```
conda activate earlgrey_env
```

Test that Earl Grey is working by bringing up the help menu

```
earlGrey -h
```

## Running Earl Grey

Run Earl Grey on the unmasked Naked Mole Rat chromosome found in `../example_data`

```
nohup earlGrey \
 -g ../example_data/NMRchr28.fa
 -s nmr_chr28 \
 -o ./example_results \
 -r eukarya \
 -d yes >& nohup.earlgrey.out
```

# Example: Generating gene models for chromosome 28 of the Naked Mole Rat

## Installations

Pull a Docker container for LiftOff and test that the container is working

```
docker pull staphb/liftoff
docker run -v "$(pwd)":/tmp staphb/liftoff liftoff -h
```

Alternatively, create a conda environment and check that liftoff is working

```
conda create -n liftoff_env liftoff=1.6.3
conda activate liftoff_env
liftoff -h
```

## Download a reference genome

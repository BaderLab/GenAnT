# How to annotate a mammalian genome

This is a tutorial on how to annotate a newly sequenced mammalian genome. This takes the user from their FASTA sequence to a high-quality GFF file annotated with gene symbols. We break the process into four main steps:

1. Repeat masking
2. Generating gene models
3. Combining and filtering gene models
4. Predicting gene function by annotating gene symbols

We recommend tools and best practices for each step, providing code to help the user execute each task. It is up to the user to install the tools that we recommend for this pipeline, however we make note of challenging installation processes that we have encountered with certain tools.

## Repeat masking

Repeat masking can be done with [Earl Grey](https://github.com/TobyBaril/EarlGrey).

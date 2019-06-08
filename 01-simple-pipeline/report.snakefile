rule report:
    message: "Generating report on the variance between samples and the reference genome."
    input:
        "calls/all.vcf"
    output:
        "report.html"
    run:
        from snakemake.utils import report

        with open(input[0]) as f:
            n_calls = sum(1 for line in f if not line.startswith("#"))

        report("""
        An example workflow
        ===================================

        Reads were mapped to the Yeas reference genome 
        and variants were called jointly with
        SAMtools/BCFtools.

        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], metadata="Author: Niels v.d. Vegt", T1=input[0])
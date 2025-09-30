@testitem "Citation Bibtex" begin
    # citation bibtex
    HEOM_buffer = IOBuffer()
    HierarchicalEOM.cite(HEOM_buffer)
    captured_output = String(take!(HEOM_buffer))

    QT_buffer = IOBuffer()
    QuantumToolbox.cite(QT_buffer)
    QuantumToolbox_output = String(take!(QT_buffer))

    @test captured_output ==
          """@article{HierarchicalEOM.jl2023,\n""" *
          """  title = {An efficient {J}ulia framework for hierarchical equations of motion in open quantum systems},\n""" *
          """  author = {Huang, Yi-Te and Kuo, Po-Chen and Lambert, Neill and Cirio, Mauro and Cross, Simon and Yang, Shen-Liang and Nori, Franco and Chen, Yueh-Nan},\n""" *
          """  journal = {Communications Physics},\n""" *
          """  publisher = {Nature Portfolio},\n""" *
          """  volume = {6},\n""" *
          """  number = {1},\n""" *
          """  pages = {313},\n""" *
          """  month = {Oct},\n""" *
          """  year = {2023},\n""" *
          """  doi = {10.1038/s42005-023-01427-2},\n""" *
          """  url = {https://doi.org/10.1038/s42005-023-01427-2}\n""" *
          """}\n""" *
          QuantumToolbox_output
end

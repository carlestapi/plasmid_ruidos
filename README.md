# plasmid_ruidos

modelo de dinámica de plásmidos en ambientes estocásticos
- cargar series temporales: { periódicas / continuas / ruidos / difusión }
- runNoise o bien runCA (para autómatas cellulares en .txt)
- cálculo: { AutoCorr / CrossCorr / CrossCov / Pearson }
- análisis frequencial: { PS / CrossCoherence }
- otros: { Entropy / JointEntropy / MutualInf / GenHurst }

generar ambientes ctes {0-1} / periódicos {sine/saw/square/envgauss} en: periodics_to_csv.py

generar ambientes noise {violet/blue/white/pink/red} / continuous stochastic / diffusion en: Conjunto_Noise_Cont_Diff.py

generar ambientes automatas celulares {cellpylib} en: continousCA.py

* crear folder /figures
* crear folder /series_temporales


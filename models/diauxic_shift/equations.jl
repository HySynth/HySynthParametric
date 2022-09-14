# constants
const k_cat1 = 3000  # max import rate of T1 (min^-1)
const k_cat2 = 2000  # max import rate of T2 (min^-1)
const KT = 1000  # enzymatic threshold (mmol)
const kr = 1260  # max elongation rate (aa min^-1)
const Kr = 7 # elongation threshold (aa min^-1)
const nQ = 300  # average length of non-ribosomal enzymes (aa)
const nRP = 300  # average length of non-ribosomal enzymes (aa)
const nT1 = 400  # length of T1 (aa)
const nT2 = 1500  # length of T2 (aa)
const nR = 7459  # length of ribosome (aa)
const w = 100  # average molar weight of amino acids (mg mmol^-1)
const kde = 0.01  # enzyme degradation rate (min^-1)
const kdRP = 0.2  # regulatory protein degradation rate (min^-1)
const kdT1 = 0.05  # degradation rate of T1 (min^-1)
const kdT2 = 0.05  # degradation rate of T2 (min^-1)

# variables notation
# x = [C1, C2, M, Q, R, T1, T2, RP]

# dynamics in the off_on mode (mode 1)
@taylorize function diauxic_off_on!(dx, x, p, t)

	# vC1 = (k_cat1 * C1 * T1)/(KT + C1)
	vC1 = (k_cat1 * x[1] * x[6])/(KT + x[1])
	# vC2 = (k_cat2 * C2 * T2)/(KT + C2)
	vC2 = (k_cat2 * x[2] * x[7])/(KT + x[2])
	# vM = (kr * M * R)/(Kr + M)
	vM = (kr * x[3] * x[5])/(Kr + x[3])
	
	# C1' = -vC1
	dx[1] = - vC1
	
	# C2' = -vC2
	dx[2] = - vC2
	
	# M' = vC1 + vC2 - vM
	dx[3] = vC1 + vC2 - vM
	
	# Q' = 1/(4 * nQ) * vM - kde * Q
	dx[4] = 1/(4 * nQ) * vM - kde * x[4]
	
	# R' = 1/(4 * nR) * vM - kde * R
	dx[5] = 1/(4 * nR) * vM - kde * x[5]
	
	# T1' = 1/(4 * nT1) * vM - kdT1 * T1
	dx[6] = 1/(4 * nT1) * vM - kdT1 * x[6]
	
	# T2' = 1/(4 * nT2) * vM - kdT2 * T2
	dx[7] = 1/(4 * nT2) * vM - kdT2 * x[7]
	
	# RP' = -kdRP * RP
	dx[8] = -kdRP * x[8]
	
	return dx
end


# dynamics in the on_on mode (mode 2)
@taylorize function diauxic_on_on!(dx, x, p, t)

	# vC1 = (k_cat1 * C1 * T1)/(KT + C1)
	vC1 = (k_cat1 * x[1] * x[6])/(KT + x[1])
	# vC2 = (k_cat2 * C2 * T2)/(KT + C2)
	vC2 = (k_cat2 * x[2] * x[7])/(KT + x[2])
	# vM = (kr * M * R)/(Kr + M)
	vM = (kr * x[3] * x[5])/(Kr + x[3])
	
	# C1' = -vC1
	dx[1] = - vC1
	
	# C2' = -vC2
	dx[2] = - vC2
	
	# M' = vC1 + vC2 - vM
	dx[3] = vC1 + vC2 - vM
	
	# Q' = 1/(5 * nQ) * vM - kde * Q
	dx[4] = 1/(5 * nQ) * vM - kde * x[4]
	
	# R' = 1/(5 * nR) * vM - kde * R
	dx[5] = 1/(5 * nR) * vM - kde * x[5]
	
	# T1' = 1/(5 * nT1) * vM - kdT1 * T1
	dx[6] = 1/(5 * nT1) * vM - kdT1 * x[6]
	
	# T2' = 1/(5 * nT2) * vM - kdT2 * T2
	dx[7] = 1/(5 * nT2) * vM - kdT2 * x[7]
	
	# RP' = 1/(5 * nRP) * vM - kdRP * RP
	dx[8] = 1/(5 * nRP) * vM - kdRP * x[8]
	
	return dx
end


# dynamics in the on_off mode (mode 3)
@taylorize function diauxic_on_off!(dx, x, p, t)

	# vC1 = (k_cat1 * C1 * T1)/(KT + C1)
	vC1 = (k_cat1 * x[1] * x[6])/(KT + x[1])
	# vC2 = (k_cat2 * C2 * T2)/(KT + C2)
	vC2 = (k_cat2 * x[2] * x[7])/(KT + x[2])
	# vM = (kr * M * R)/(Kr + M)
	vM = (kr * x[3] * x[5])/(Kr + x[3])
	
	# C1' = -vC1
	dx[1] = - vC1
	
	# C2' = -vC2
	dx[2] = - vC2
	
	# M' = vC1 + vC2 - vM
	dx[3] = vC1 + vC2 - vM
	
	# Q' = 1/(4 * nQ) * vM - kde * Q
	dx[4] = 1/(4 * nQ) * vM - kde * x[4]
	
	# R' = 1/(4 * nR) * vM - kde * R
	dx[5] = 1/(4 * nR) * vM - kde * x[5]
	
	# T1' = 1/(4 * nT1) * vM - kdT1 * T1
	dx[6] = 1/(4 * nT1) * vM - kdT1 * x[6]
	
	# T2' = - kdT2 * T2
	dx[7] = - kdT2 * x[7]
	
	# RP' = 1/(4 * nRP) * vM - kdRP * RP
	dx[8] = 1/(4 * nRP) * vM - kdRP * x[8]
	
	return dx
end


# dynamics in the off_off mode (mode 4)
@taylorize function diauxic_off_off!(dx, x, p, t)

	# vC1 = (k_cat1 * C1 * T1)/(KT + C1)
	vC1 = (k_cat1 * x[1] * x[6])/(KT + x[1])
	# vC2 = (k_cat2 * C2 * T2)/(KT + C2)
	vC2 = (k_cat2 * x[2] * x[7])/(KT + x[2])
	# vM = (kr * M * R)/(Kr + M)
	vM = (kr * x[3] * x[5])/(Kr + x[3])
	
	# C1' = -vC1
	dx[1] = - vC1
	
	# C2' = -vC2
	dx[2] = - vC2
	
	# M' = vC1 + vC2 - vM
	dx[3] = vC1 + vC2 - vM
	
	# Q' = 1/(3 * nQ) * vM - kde * Q
	dx[4] = 1/(3 * nQ) * vM - kde * x[4]
	
	# R' = 1/(3 * nR) * vM - kde * R
	dx[5] = 1/(3 * nR) * vM - kde * x[5]
	
	# T1' = 1/(3 * nT1) * vM - kdT1 * T1
	dx[6] = 1/(3 * nT1) * vM - kdT1 * x[6]
	
	# T2' = - kdT2 * T2
	dx[7] = - kdT2 * x[7]
	
	# RP' = - kdRP * RP
	dx[8] = - kdRP * x[8]
	
	return dx
end


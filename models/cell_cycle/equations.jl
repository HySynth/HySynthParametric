# constants
const ksa1 = 5  	# parameter for synthesis
const ksa2 = 6		# parameter for synthesis
const ksa3 = 20		# parameter for synthesis
const kda1 = 0.2	# parameter for degradation
const kda2 = 1.2	# parameter for degradation
const kda3 = 1.2	# parameter for degradation
const ksb1 = 2.5	# parameter for synthesis
const ksb2 = 6		# parameter for synthesis
const kdb1 = 0.2	# parameter for degradation
const kdb2 = 1.2	# parameter for degradation
const kdb3 = 0.3	# parameter for degradation
const kse1 = 0.02	# parameter for synthesis
const kse2 = 2		# parameter for synthesis
const kde1 = 0.02	# parameter for degradation (misprint in Table 1?)
const kde2 = 0.5 	# parameter for degradation
const γ = 0.029		# specific growth rate

# variables notation
# x = [CyCA, CyCB, CyCE, M, t]

# dynamics in the G1a mode (mode 1)
function mode_G1a()
	# boolean variables
	BTFE = 0	# transcription factor TFE
	BSCF = 0	# ubiquitin-ligase SCF
	BTFB = 0	# transcription factor TFB
	BCdc20A = 0	# Cdc20 active on cyclin A
	BCdc20B = 0	# Cdc20 active on cyclin B
	BCdh1 = 1	# Cdh1 active
	
	# rate for synthesis of cyclin A
	ksa = ksa1 + ksa2 * BTFE + ksa3 * BTFB
	# rate for degradation of cyclin A
	kda = kda1 + kda2 * BCdc20A + kda3 * BCdh1
	# rate for synthesis of cyclin B
	ksb = ksb1 + ksb2 * BTFB
	# rate for degradation of cyclin B
	kdb = kdb1 + kdb2 * BCdc20B + kdb3 * BCdh1
	# rate for synthesis of cyclin E
	kse = kse1 + kse2 * BTFE
	# rate for degradation of cyclin E
	kde = kde1 + kde2 * BSCF
	
	A = zeros(5, 5)
	b = zeros(5)
	
	# CyCA' = ksa - kda * CyCA
	A[1, 1] = -kda
	b[1] = ksa
	
	# CyCB' = ksb - kdb * CyCB
	A[2, 2] = -kdb
	b[2] = ksb
	
	# CyCE' = kse - kde * CyCE
	A[3, 3] = -kde
	b[3] = kse
	
	# M' = γ * M
	A[4, 4] = γ
	
	# t' = 1
	b[5] = 1
	
	return A, b
end

# dynamics in the early_G1b mode (mode 2)
function mode_early_G1b()

	# boolean variables
	BTFE = 1	# transcription factor TFE
	BSCF = 0	# ubiquitin-ligase SCF
	BTFB = 0	# transcription factor TFB
	BCdc20A = 0	# Cdc20 active on cyclin A
	BCdc20B = 0	# Cdc20 active on cyclin B
	BCdh1 = 1	# Cdh1 active
	
	# rate for synthesis of cyclin A
	ksa = ksa1 + ksa2 * BTFE + ksa3 * BTFB
	# rate for degradation of cyclin A
	kda = kda1 + kda2 * BCdc20A + kda3 * BCdh1
	# rate for synthesis of cyclin B
	ksb = ksb1 + ksb2 * BTFB
	# rate for degradation of cyclin B
	kdb = kdb1 + kdb2 * BCdc20B + kdb3 * BCdh1
	# rate for synthesis of cyclin E
	kse = kse1 + kse2 * BTFE
	# rate for degradation of cyclin E
	kde = kde1 + kde2 * BSCF
	
	A = zeros(5, 5)
	b = zeros(5)
	
	# CyCA' = ksa - kda * CyCA
	A[1, 1] = -kda
	b[1] = ksa
	
	# CyCB' = ksb - kdb * CyCB
	A[2, 2] = -kdb
	b[2] = ksb
	
	# CyCE' = kse - kde * CyCE
	A[3, 3] = -kde
	b[3] = kse
	
	# M' = γ * M
	A[4, 4] = γ
	
	# t' = 1
	b[5] = 1
	
	return A, b
end


# dynamics in the late_G1b mode (mode 3)
function mode_late_G1b()

	# boolean variables
	BTFE = 1	# transcription factor TFE
	BSCF = 0	# ubiquitin-ligase SCF
	BTFB = 0	# transcription factor TFB
	BCdc20A = 0	# Cdc20 active on cyclin A
	BCdc20B = 0	# Cdc20 active on cyclin B
	BCdh1 = 0	# Cdh1 active
	
	# rate for synthesis of cyclin A
	ksa = ksa1 + ksa2 * BTFE + ksa3 * BTFB
	# rate for degradation of cyclin A
	kda = kda1 + kda2 * BCdc20A + kda3 * BCdh1
	# rate for synthesis of cyclin B
	ksb = ksb1 + ksb2 * BTFB
	# rate for degradation of cyclin B
	kdb = kdb1 + kdb2 * BCdc20B + kdb3 * BCdh1
	# rate for synthesis of cyclin E
	kse = kse1 + kse2 * BTFE
	# rate for degradation of cyclin E
	kde = kde1 + kde2 * BSCF
	
	A = zeros(5, 5)
	b = zeros(5)
	
	# CyCA' = ksa - kda * CyCA
	A[1, 1] = -kda
	b[1] = ksa
	
	# CyCB' = ksb - kdb * CyCB
	A[2, 2] = -kdb
	b[2] = ksb
	
	# CyCE' = kse - kde * CyCE
	A[3, 3] = -kde
	b[3] = kse
	
	# M' = γ * M
	A[4, 4] = γ
	
	# t' = 1
	b[5] = 1
	
	return A, b
end


# dynamics in the S mode (mode 4)
function mode_S()

	# boolean variables
	BTFE = 1	# transcription factor TFE
	BSCF = 1	# ubiquitin-ligase SCF
	BTFB = 0	# transcription factor TFB
	BCdc20A = 0	# Cdc20 active on cyclin A
	BCdc20B = 0	# Cdc20 active on cyclin B
	BCdh1 = 0	# Cdh1 active
	
	# rate for synthesis of cyclin A
	ksa = ksa1 + ksa2 * BTFE + ksa3 * BTFB
	# rate for degradation of cyclin A
	kda = kda1 + kda2 * BCdc20A + kda3 * BCdh1
	# rate for synthesis of cyclin B
	ksb = ksb1 + ksb2 * BTFB
	# rate for degradation of cyclin B
	kdb = kdb1 + kdb2 * BCdc20B + kdb3 * BCdh1
	# rate for synthesis of cyclin E
	kse = kse1 + kse2 * BTFE
	# rate for degradation of cyclin E
	kde = kde1 + kde2 * BSCF
	
	A = zeros(5, 5)
	b = zeros(5)
	
	# CyCA' = ksa - kda * CyCA
	A[1, 1] = -kda
	b[1] = ksa
	
	# CyCB' = ksb - kdb * CyCB
	A[2, 2] = -kdb
	b[2] = ksb
	
	# CyCE' = kse - kde * CyCE
	A[3, 3] = -kde
	b[3] = kse
	
	# M' = γ * M
	A[4, 4] = γ
	
	# t' = 1
	b[5] = 1
	
	return A, b
end


# dynamics in the G2 mode (mode 5)
function mode_G2()

	# boolean variables
	BTFE = 1	# transcription factor TFE
	BSCF = 1	# ubiquitin-ligase SCF
	BTFB = 1	# transcription factor TFB
	BCdc20A = 0	# Cdc20 active on cyclin A
	BCdc20B = 0	# Cdc20 active on cyclin B
	BCdh1 = 0	# Cdh1 active
	
	# rate for synthesis of cyclin A
	ksa = ksa1 + ksa2 * BTFE + ksa3 * BTFB
	# rate for degradation of cyclin A
	kda = kda1 + kda2 * BCdc20A + kda3 * BCdh1
	# rate for synthesis of cyclin B
	ksb = ksb1 + ksb2 * BTFB
	# rate for degradation of cyclin B
	kdb = kdb1 + kdb2 * BCdc20B + kdb3 * BCdh1
	# rate for synthesis of cyclin E
	kse = kse1 + kse2 * BTFE
	# rate for degradation of cyclin E
	kde = kde1 + kde2 * BSCF
	
	A = zeros(5, 5)
	b = zeros(5)
	
	# CyCA' = ksa - kda * CyCA
	A[1, 1] = -kda
	b[1] = ksa
	
	# CyCB' = ksb - kdb * CyCB
	A[2, 2] = -kdb
	b[2] = ksb
	
	# CyCE' = kse - kde * CyCE
	A[3, 3] = -kde
	b[3] = kse
	
	# M' = γ * M
	A[4, 4] = γ
	
	# t' = 1
	b[5] = 1
	
	return A, b
end


# dynamics in the prophase mode (mode 6)
function mode_prophase()

	# boolean variables
	BTFE = 0	# transcription factor TFE
	BSCF = 1	# ubiquitin-ligase SCF
	BTFB = 1	# transcription factor TFB
	BCdc20A = 0	# Cdc20 active on cyclin A
	BCdc20B = 0	# Cdc20 active on cyclin B
	BCdh1 = 0	# Cdh1 active
	
	# rate for synthesis of cyclin A
	ksa = ksa1 + ksa2 * BTFE + ksa3 * BTFB
	# rate for degradation of cyclin A
	kda = kda1 + kda2 * BCdc20A + kda3 * BCdh1
	# rate for synthesis of cyclin B
	ksb = ksb1 + ksb2 * BTFB
	# rate for degradation of cyclin B
	kdb = kdb1 + kdb2 * BCdc20B + kdb3 * BCdh1
	# rate for synthesis of cyclin E
	kse = kse1 + kse2 * BTFE
	# rate for degradation of cyclin E
	kde = kde1 + kde2 * BSCF
	
	A = zeros(5, 5)
	b = zeros(5)
	
	# CyCA' = ksa - kda * CyCA
	A[1, 1] = -kda
	b[1] = ksa
	
	# CyCB' = ksb - kdb * CyCB
	A[2, 2] = -kdb
	b[2] = ksb
	
	# CyCE' = kse - kde * CyCE
	A[3, 3] = -kde
	b[3] = kse
	
	# M' = γ * M
	A[4, 4] = γ
	
	# t' = 1
	b[5] = 1
	
	return A, b
end


# dynamics in the metaphase mode (mode 7)
function mode_metaphase()

	# boolean variables
	BTFE = 0	# transcription factor TFE
	BSCF = 1	# ubiquitin-ligase SCF
	BTFB = 1	# transcription factor TFB
	BCdc20A = 1	# Cdc20 active on cyclin A
	BCdc20B = 0	# Cdc20 active on cyclin B
	BCdh1 = 0	# Cdh1 active
	
	# rate for synthesis of cyclin A
	ksa = ksa1 + ksa2 * BTFE + ksa3 * BTFB
	# rate for degradation of cyclin A
	kda = kda1 + kda2 * BCdc20A + kda3 * BCdh1
	# rate for synthesis of cyclin B
	ksb = ksb1 + ksb2 * BTFB
	# rate for degradation of cyclin B
	kdb = kdb1 + kdb2 * BCdc20B + kdb3 * BCdh1
	# rate for synthesis of cyclin E
	kse = kse1 + kse2 * BTFE
	# rate for degradation of cyclin E
	kde = kde1 + kde2 * BSCF
	
	A = zeros(5, 5)
	b = zeros(5)
	
	# CyCA' = ksa - kda * CyCA
	A[1, 1] = -kda
	b[1] = ksa
	
	# CyCB' = ksb - kdb * CyCB
	A[2, 2] = -kdb
	b[2] = ksb
	
	# CyCE' = kse - kde * CyCE
	A[3, 3] = -kde
	b[3] = kse
	
	# M' = γ * M
	A[4, 4] = γ
	
	# t' = 1
	b[5] = 1
	
	return A, b
end


# dynamics in the anaphase mode (mode 8)
function mode_anaphase()

	# boolean variables
	BTFE = 0	# transcription factor TFE
	BSCF = 1	# ubiquitin-ligase SCF
	BTFB = 1	# transcription factor TFB
	BCdc20A = 1	# Cdc20 active on cyclin A
	BCdc20B = 1	# Cdc20 active on cyclin B
	BCdh1 = 0	# Cdh1 active
	
	# rate for synthesis of cyclin A
	ksa = ksa1 + ksa2 * BTFE + ksa3 * BTFB
	# rate for degradation of cyclin A
	kda = kda1 + kda2 * BCdc20A + kda3 * BCdh1
	# rate for synthesis of cyclin B
	ksb = ksb1 + ksb2 * BTFB
	# rate for degradation of cyclin B
	kdb = kdb1 + kdb2 * BCdc20B + kdb3 * BCdh1
	# rate for synthesis of cyclin E
	kse = kse1 + kse2 * BTFE
	# rate for degradation of cyclin E
	kde = kde1 + kde2 * BSCF
	
	A = zeros(5, 5)
	b = zeros(5)
	
	# CyCA' = ksa - kda * CyCA
	A[1, 1] = -kda
	b[1] = ksa
	
	# CyCB' = ksb - kdb * CyCB
	A[2, 2] = -kdb
	b[2] = ksb
	
	# CyCE' = kse - kde * CyCE
	A[3, 3] = -kde
	b[3] = kse
	
	# M' = γ * M
	A[4, 4] = γ
	
	# t' = 1
	b[5] = 1
	
	return A, b
end


# dynamics in the telophase mode (mode 9)
function mode_telophase()

	# boolean variables
	BTFE = 0	# transcription factor TFE
	BSCF = 1	# ubiquitin-ligase SCF
	BTFB = 0	# transcription factor TFB
	BCdc20A = 1	# Cdc20 active on cyclin A
	BCdc20B = 1	# Cdc20 active on cyclin B
	BCdh1 = 1	# Cdh1 active
	
	# rate for synthesis of cyclin A
	ksa = ksa1 + ksa2 * BTFE + ksa3 * BTFB
	# rate for degradation of cyclin A
	kda = kda1 + kda2 * BCdc20A + kda3 * BCdh1
	# rate for synthesis of cyclin B
	ksb = ksb1 + ksb2 * BTFB
	# rate for degradation of cyclin B
	kdb = kdb1 + kdb2 * BCdc20B + kdb3 * BCdh1
	# rate for synthesis of cyclin E
	kse = kse1 + kse2 * BTFE
	# rate for degradation of cyclin E
	kde = kde1 + kde2 * BSCF
	
	A = zeros(5, 5)
	b = zeros(5)
	
	# CyCA' = ksa - kda * CyCA
	A[1, 1] = -kda
	b[1] = ksa
	
	# CyCB' = ksb - kdb * CyCB
	A[2, 2] = -kdb
	b[2] = ksb
	
	# CyCE' = kse - kde * CyCE
	A[3, 3] = -kde
	b[3] = kse
	
	# M' = γ * M
	A[4, 4] = γ
	
	# t' = 1
	b[5] = 1
	
	return A, b
end

G1a! = mode_G1a()
early_G1b! = mode_early_G1b()
late_G1b! = mode_late_G1b()
S! = mode_S()
G2! = mode_G2()
prophase! = mode_prophase()
metaphase! = mode_metaphase()
anaphase! = mode_anaphase()
telophase! = mode_telophase()

HOW TO FILL IN EXCEL FILE


Important remarks:
- If your experiment includes the hybridization of two libraries, you have to fill in two separate Excel files (i.e., one file per library).
- Do only edit the sheets "step1", ..., "stepN", "scaffolds", and "const" (i.e., the sheet "smarts" can be ignored).
- Delete the sheets you don't need (e.g., delete "step2" if your experiment includes only one reaction step).
- Do only use the options available (see below).
- In case of a two step reaction (e.g., SR followed by ABF), write all reaction types separated by a comma (see example usage).
- If no reaction takes place or no scaffold is used, leave the corresponding cell empty (see example usage).
- Only reactions between two building blocks or between a building block and a scaffold are considered (i.e., reactions between a building block and a DNA strand should be neglected).
- All sequences provided should be written in the 5'-to-3' direction.
- In case of a hybridization, the constant regions provided must not overlap (see example usage).



##################
# STEP 1, ..., N #
##################

ID
Description:	ID of the N-th building block
Options:	1, ..., n

SMILES
Description:	SMILES of building block n
Options:	Valid SMILES structure

ScaffoldID
Description:	ID of the corresponding scaffold (only for step 1)
Options		1, ..., m

Codon
Description:	Codon of building block n
Option:		[ATCG]

ReactionType
Description:	Reaction type(s) (amide bond formation, cycloaddition, dehalogenation, diazo transfer, intermolecular nucleophilic substitution (type 2), intramolecular nucleophilic substitution (type 2), Sonogashira, Staudinger reduction, Suzuki)
Options:	ABF, CuAAC, DH, DT, SN2-2, SN2-1, Son, SR, Suz


Example usage:

ID	SMILES		ScaffoldID	Codon	ReactionType
------------------------------------------------------------
1	NCC#C		1		GTGGTG	CuAAC
2	OCC(O)=O	1		ATATTG	SR, ABF
3	[I]		2		TGGCAA



#############
# SCAFFOLDS #
#############

ScaffoldID
Description:	ID of the corresponding scaffold
Options:	1, ..., m

SMILES
Description:	SMILES of scaffold m
Options:	Valid SMILES structure


Example usage:

ScaffoldID	SMILES
------------------------------------------------
1		Ic1ccc(CC(N=[N+]=[N-])C(O)=O)cc1
2		[N-]=[N+]=NC(C(O)=O)Cc1cc(I)ccc1



#########
# CONST #
#########

Sequence
Description:	Sequence of the constant regions surrounding the codons of the building blocks
Options:	[ATCG]

Reverse
Description:	Whether the final DNA sequence contains the reverse sequence of the constant region provided
Option:		0, 1

Complement
Description:	Whether the final DNA sequence contains the complement sequence of the constant region provided
Options:	0, 1


Example usage (one library):

Final DNA sequence:	5' GAC{codon}GCT{codon}GCT 3'

Sequence			Reverse	Complement
--------------------------------------------------
GAC{codon}GCT{codon}GCT		0	0


Example usage (hybridization of two libraries):

5' GTA{codon}GGA              3'
   |||       |||
3' CAT{xxxxx}CCTGAA{codon}TGC 5'

Final DNA sequence:	5' GTA{codon}GGACTT{codon}ACG 3'

Excel file #1:
Sequence				Reverse	Complement
----------------------------------------------------------
GTA{codon}GGA				0	0

Excel file #2:
Sequence				Reverse	Complement
----------------------------------------------------------
CGT{codon}AAG				1	1


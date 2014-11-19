//-------|---------|---------|---------|---------|---------|---------|---------|
/*
file_format.h - 
*/

#ifndef	FILE_FORMAT_H
#define	FILE_FORMAT_H

#include "all.h"

/*
COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angst
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angst
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angst
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
*/
typedef struct
{
	char		RecName[6];		// 1  - 6
	char		UID[5];			// 7  - 11
	char		Blank;			// 12
	char		Name[4];		// 13 - 16
	char		AltLoc;			// 17
	char		ResName[3];		// 18 - 20
	char		Blank2;			// 21
	char		ChainId;		// 22
	char		ResSeq[4];		// 23 - 26
	char		ICode;			// 27
	char		Blank3[3];		// 28 - 30
	char		X[8];			// 31 - 38
	char		Y[8];			// 39 - 46
	char		Z[8];			// 47 - 54
	char		Occupancy[6];	// 55 - 60
	char		TempFactor[6];	// 61 - 66
	char		Blank4[10];		// 67 - 76
	char		Element[2];		// 77 - 78
	char		Charge[2];		// 79 - 80
} FILE_FORMAT_PDB_RECORD;

typedef union
{
	char					Buffer[80+1];
	FILE_FORMAT_PDB_RECORD	Atom;
} FILE_FORMAT_PDB;

typedef struct
{
	char				NumRecs[8];		//	1  - 8
	char				Blank1;			//	9
	char				Type[6];		//	10 - 15
} FILE_FORMAT_PSF_HEADER;

typedef struct	//	NOTE: this format is essentially a guess.
{
	char				UID[8];			//	1  - 8
	char				Blank1;			//	9
	char				Unk1[4];		//	10 - 13
	char				Blank2;			//	14
	char				Unk2[4];		//	15 - 18
	char				Blank3;			//	19
	char				Unk3[4];		//	20 - 23
	char				Blank4;			//	24
	char				Unk4[4];		//	25 - 28
	char				Blank5;			//	29
	char				Unk5[4];		//	30 - 33
	char				Blank6;			//	34
	char				Charge[16];		//	35 - 50
	char				Blank7;			//	51
	char				Mass[7];		//	52 - 58
	char				Blank8;			//	59
	char				Unk6[11];		//	60 - 70
} FILE_FORMAT_PSF_RECORD;

typedef union
{
	char					Buffer[80+1];
	FILE_FORMAT_PSF_HEADER	Header;
	FILE_FORMAT_PSF_RECORD	Atom;
} FILE_FORMAT_PSF;

#endif

// End of file
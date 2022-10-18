# Deduper: 
Identify and remove PCR duplicates from a SAM file
---
**SAM FILE Example record:**\
(Input: Sorted SAM file | Output: Unique entry SAM file (no pcr duplicates))
```
NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	2	76814284	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
```
**PCR Duplicates:** 
- Same alignment position (Chromosome, position, strand)
- Potential Soft clipping
- Same Unique Molecular Index (UMI)
- We will have single-ended reads for this assignment

**Fields of interest:**
- RNAME (SAM col 3)
- POS (SAM col 4)
- FLAG (SAM col 2)
- CIGAR (SAM col 6)
- QNAME "UMI" (SAM col 1) 

**Concepts:**
- Bitflag: 16 of FLAG will tell us if it's strand specific.\
![image](https://user-images.githubusercontent.com/106117735/195733004-b9bc7d2b-81e1-4218-a6bc-ec6c65170863.png)

- Soft Clipping: occurs at the ends of an alignment. Noted for in CIGAR string.\
![image](https://user-images.githubusercontent.com/106117735/195733139-99a02685-c427-46a5-ba3a-745d861d1703.png)


**Functions:**\
  def recorder(line) -> lst
  
  "Reads a seq record and pulls fields of interest as set" 
  input: line
  output: (RNAME,POS,FLAG, CIGAR, QNAME)
  ```
  input : regex the line and pull RNAME, POS, QNAME...etc.
  output: set of elements pulled from record
  ```
  def strandflag(flag: int) -> str:
  
  "Decodes bitwise flag, returns str of "1" or "0" to check for strand specificity."
   input: 36
   output: 1
   ```
   binary = bin(flag) 
   strand_spec = str(binary)[-5] 
   return(strand_spec)
   ```
  def clipcig(cig: str) -> int:
    
  "Interprets cigar string, and returns the sum of soft clipped(S) and matches (M) May need to apply to POS? (should be leftmost)"\
  input: 71M\
  output: 71\
  input: 3S68M\
  output 71
  
  NOTE: For position correction, take scliplist sum and subtract from POS.\
  input: pos = 234 (given 3S clipped bases)\
  output: corrected = 231 
  ```
  scliplist = []
  matchlist = []
 append scliplist int re.findall(\d+)S
 append matchlist int re.findall(\d+)M
 abslen = sum(scliplist) + sum(matchlist) 
 return(abslen)
 ```

 
 
 **Order of Operations:**
 - SORT Data by chromosome?
 - Store known UMIs
 - Loop through lines in file. (Applying conditional funnel and removal steps) 
 - Conditional funnel: (most broad to specific) *using recorder() result*
    1. IF chromosome name (RNAME) is the same
    2. IF strand is the same *use strandflag()* 
    3. IF position is the same given potential clipping *use clipcig()*
    4. IF UMI is known and is the same *could potentially write function?*
 - IF UMI is unknown, interpret as *erroneous* given known UMIs are used, store these elsewhere? 
 - If funnel works and runs to completion, surviving seqs should be PCR duplicates. 
 - Omit the duplicates, and store elsewhere? 
 - For performance, may chunk the workload by chromosome?
 
 Note: Negative strand adjusted position must be calculated differently.
 Will output to a new deduped SAM file.

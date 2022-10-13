# Deduper: 
Identify and remove PCR duplicates from a SAM file. 
---
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
- Bitflag: 16 of FLAG will tell us if it's strand specific. 
- Soft Clipping: occurs at the ends of an alignment. Noted for in CIGAR string.

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
    
  "Interprets cigar string, and returns the sum of soft clipped(S) and matches (M) May need to apply to POS? (should be leftmost)"
  input: 71M
  output: 71
  input: 3S68M
  output 71
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
 - Conditional funnel: (most broad to specific) *using recorder() result*
    1. IF chromosome name (RNAME) is the same
    2. IF strand is the same *use strandflag()* 
    3. IF position is the same given potential clipping *use clipcig()*
    4. IF UMI is known and is the same *could potentially write function?*
 - IF UMI is unknown, interpret as *erroneous* given known UMIs are used, store these elsewhere? 
 - If funnel works and runs to completion, surviving seqs should be PCR duplicates. 
 - Omit the duplicates, and store elsewhere? 
 
 
  

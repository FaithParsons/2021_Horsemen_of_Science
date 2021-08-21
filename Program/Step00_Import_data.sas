

proc import datafile= "..\data\HR BEHAV database for Faith.xlsx"
	/* only read rows where timepoint is 2, 3, or 4*/
	out = mydata (WHERE=(timepoint IN (2, 3, 4))) 
	dbms = xlsx replace;
	sheet = "data";
run;

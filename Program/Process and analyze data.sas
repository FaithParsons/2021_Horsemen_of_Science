* Import data;
proc import datafile= "..\data\HR BEHAV database for Faith.xlsx"
	/* only read rows where timepoint is 2, 3, or 4*/
	out = import (WHERE=(timepoint IN (2, 3, 4))) 
	dbms = xlsx replace;
	sheet = "data";
run;


data horses;
	set import;
	* Create a treatment variable that identifies the type of treatment;
	length treatment $5;
	if group = 1 and timepoint = 3 then treatment = "cow";
	else if group = 2 and timepoint = 3 then treatment = "lion";
	else treatment = "food"; 

	* Create a group variable;
	if group = 1 THEN group = "A";
	else if group = 2 THEN group = "B";
	
	* Correct the data entry issue;
	if approach_lat > Visible then approach_lat_corrected = Visible;
	else approach_lat_corrected  = approach_lat;

	pct_approach_lat = 100*approach_lat_corrected /visible;

run;

title1 "Outcome: sniff_dur/visible";
proc glimmix data=mydata;
class treatment (ref='cow') id timepoint (ref='2') group (ref='A');
model sniff_dur/visible = treatment     
	/ solution dist=binomial oddsratio;
random intercept group pct_latency / subject=id  ; 
run; 

title1 "Outcome: Ears_forw/visible";
proc glimmix data=mydata;
class treatment (ref='cow') id timepoint (ref='2') group (ref='A');
model Ears_forw/visible = treatment     
	/ solution dist=binomial oddsratio;
random intercept group pct_latency / subject=id  ; 
run; 


title1 "Outcome: withdr_dur/visible";
proc glimmix data=mydata;
class treatment (ref='cow') id timepoint (ref='2') group (ref='A');
model withdr_dur/visible = treatment     
	/ solution dist=binomial oddsratio;
random intercept group pct_latency / subject=id  ; 
run; 


title1 "Outcome: Interact_dur/visible";
proc glimmix data=mydata;
class treatment (ref='cow') id group (ref='A');
model Interact_dur/visible = treatment     
	/ solution dist=binomial oddsratio;
random  _residual_ / subject=id; 
run; 


title1 "Outcome: Interact_dur/visible";
proc glimmix data=mydata;
class treatment (ref='cow') id group (ref='A');
model Interact_dur/visible = treatment     
	/ solution dist=binomial oddsratio;
random  group / subject=id type=ar(1); 
run; 


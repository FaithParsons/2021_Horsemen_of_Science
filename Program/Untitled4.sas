
proc glimmix data=mydata;
class treatment (ref='cow') id timepoint (ref='2') stimulus (ref='1');
model sniff_dur/visible = treatment     
	/ solution dist=binomial;
random intercept stimulus pct_latency / subject=id  ; 
run; 


proc glimmix data=mydata;
class treatment (ref='cow') id timepoint (ref='2') stimulus (ref='1');
model Ears_forw/visible = treatment     
	/ solution dist=binomial;
random intercept stimulus pct_latency / subject=id  ; 
run; 


proc glimmix data=mydata;
class treatment (ref='cow') id timepoint (ref='2') stimulus (ref='1');
model withdr_dur/visible = treatment     
	/ solution dist=binomial;
random intercept stimulus  pct_latency/ subject=id  ; 
run; 

proc glimmix data=mydata;
class treatment (ref='cow') id timepoint (ref='2') stimulus (ref='1');
model Interact_dur/visible = treatment     
	/ solution dist=binomial;
random _residual_ / subject=id  ; 
run; 


proc glimmix data=mydata;
class treatment (ref='cow') id timepoint (ref='2') stimulus (ref='1');
model blows_occur = treatment     
	/ solution dist=poisson;
random _residual_/ subject=id  ; 
run; 



proc sql;
select distinct  stimulus, breed, count(distinct id) as cntIDs
from mydata 
where sniff_dur ne .
group by 1, 2;
quit;


proc sql;
select distinct breed, count(distinct id) as cntIDs
from mydata 
where sniff_dur ne .
group by 1;
quit;



proc sql;
select distinct age, count(distinct id) as cntIDs
from mydata 
where sniff_dur ne .
group by 1;
quit;

proc sgpanel data=mydata(where=(timepoint in (2,3,4)));
panelby timepoint;
vbox meanhr/category=stimulus;
run; 

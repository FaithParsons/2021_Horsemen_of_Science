
* Import data;
proc import datafile= "..\data\HR BEHAV database for Faith.xlsx"
	/* only read rows where timepoint is 2, 3, or 4*/
	out = import (WHERE=(timepoint IN (2, 3, 4))) 
	dbms = xlsx replace;
	sheet = "data";
run;
option mprint;
/*proc contents data=mydata varnum; run;*/

data mydata;
	set import;
	pct_sniff = sniffpc/100;
	pct_chew = chewpc/100;
	pct_ears = earspc/100;
	pct_withdraw = WithdrPC/100;
	pct_interact = interactpc/100;
	length treatment $5;
	if stimulus = 1 and timepoint = 3 then treatment = "cow";
	else if stimulus = 2 and timepoint = 3 then treatment = "lion";
	else treatment = "food"; 
	if stimulus = 1 THEN group = "A";
	else if stimulus = 2 THEN group = "B";
	
	if stimulus = 1 and timepoint = 3 then cow = 1; else cow=0;
	if stimulus = 2 and timepoint = 3 then lion = 1; else lion = 0;
	if timepoint in (2, 4) then food = 1; else food=0;

	if approach_lat > Visible then latency_corrected = Visible;
	else latency_corrected = approach_lat;

	pct_latency = 100*latency_corrected /visible;

run;

proc freq data=mydata;
	table timepoint*stimulus*cow*lion*food/list missing;
run;

*BETA REGReSSION;
%Macro Beta_Regression(
	Dataset,
	tech,
	details,
	mu_vars,
	phi_vars,
	zero_vars,
	one_vars,depvar);
%if &zero_vars ne and &one_vars ne %then 
	%Beta_Regression_Zero_One(&Dataset,&tech,&details,&mu_vars,&phi_vars,&zero_vars,&one_vars,&depvar);
%else %if &zero_vars ne %then 
	%Beta_Regression_Zero(&Dataset,&tech,&details,&mu_vars,&phi_vars,&zero_vars,&depvar);
%else %if &one_vars ne %then 
	%Beta_Regression_One(&Dataset,&tech,&details,&mu_vars,&phi_vars,&one_vars,&depvar);
/*%else */
/*	%Beta_Regression_Only(&Dataset,&tech,&details,&mu_vars,&phi_vars,&depvar);*/
%mend;


%Macro Preprocessing(Vars,b0,xb2,b); 
data HPG;    
	%global &xb2;    
	length  &xb2 $200.;    
	&xb2=&b0;

	%if &Vars ne '' %then %do; 
		%let n=1;          
		var&n="%scan(&Vars,&n,' ')"; 
		%do %while( %scan(&Vars,&n,' ') ne  ); 
			%let n=%eval(&n+1);                   
			var&n="%scan(&Vars,&n,' ')"; 
		%end; 
		%let n_1=%eval(&n-1);          

		array xbv {*} $ var1--var&n_1;  

		%do j=1 %to &n_1;                   
			&b&j=      "&b&j"; 
		%end; %let one=1;          
		array &b{*} $8 &b&one--&b&n_1 ;           
		array   p{1} $ 8 ('+');          
		array   m{1} $ 8 ('*');          

		do   i=1 to dim(xbv) while (xbv{i} ne '');                
			&xb2= cats(of &xb2 p{1} &b{i} m{1} xbv{i});          
		end;      
	%end;    
	call   symput("&xb2",&xb2); 
run; 
%mend;

%Macro Postprocessing(hat,predict,depvar);
data &predict;
	retain record &depvar &hat;
	set &predict;
	if &depvar = . then &hat = .;
		else &hat=pred;
	record=_n_;
	keep record &depvar &hat;
run;
%mend;

%Macro Beta_Regression_Zero(
	Dataset,
	tech,
	details,
	mu_vars,
	phi_vars,
	zero_vars,
	depvar);

	%Preprocessing(&mu_vars, 'b0'  ,xb,b);
	%Preprocessing(&phi_vars,'d0'  ,wd,d);
	%Preprocessing(&zero_vars,'zero0',zeroxb,zero);

	proc nlmixed data= &Dataset tech=&tech &details;
		pizero = exp(&zeroxb)/(1 + exp(&zeroxb));
		mu=exp(&xb)/(1+exp(&xb));
		phi=exp(&wd);
		w=mu*phi;
		t=phi - mu*phi;
		if(&depvar=0) then ll=log(pizero);
			else ll=lgamma(w+t) - lgamma(w) - lgamma(t)+((w-1)*log(&depvar)) + 
				((t-1)*log(1 - &depvar))+log(1-pizero);
		model &depvar ~ general(ll);
		predict mu out=mu_results (keep = &depvar pred);
		predict phi out=phi_results(keep=&depvar pred);
		predict pizero out=pizero_results(keep=&depvar pred);
	run;
	%Postprocessing(mu_hat,mu_results,&depvar);
	%Postprocessing(phi_hat,phi_results,&depvar);
	%Postprocessing(pizero_hat,pizero_results,&depvar);

	data prediction;
		merge mu_results phi_results pizero_results;
		by record;
	run;
%mend;


%Macro Beta_Regression_One(
	Dataset,
	tech,
	details,
	mu_vars,
	phi_vars,
	one_vars,
	depvar);

	%Preprocessing(&mu_vars,'b0'  ,xb,b);
	%Preprocessing(&phi_vars,'d0'  ,wd,d);
	%Preprocessing(&one_vars,'one0',onexb,one);
	proc nlmixed data=&Dataset tech=&tech &details;
		pione=exp(&onexb)/(1+exp(&onexb));
		mu=exp(&xb)/(1+exp(&xb));
		phi=exp(&wd);
		w=mu*phi;
		t=phi - mu*phi;
		if(&depvar=1) then ll=log(pione);
			else ll=lgamma(w+t) - lgamma(w) - lgamma(t)+((w-1)*log(&depvar)) + 
				((t-1)*log(1 - &depvar)) + log(1 - pione);
		model &depvar ~ general(ll);
		predict mu out=mu_results(keep=&depvar pred);
		predict phi out=phi_results(keep=&depvar pred);
		predict pione out=pione_results(keep=&depvar pred);
	run;
	%Postprocessing(mu_hat,mu_results,&depvar);
	%Postprocessing(phi_hat,phi_results,&depvar);
	%Postprocessing(pione_hat,pione_results,&depvar);

	data prediction;
		merge mu_results phi_results pione_results;
 		by record;
	run;
%mend;

%Macro Beta_Regression_Zero_One(
	Dataset,
	tech,
	details,
	mu_vars,
	phi_vars,
	zero_vars,
	one_vars,
	depvar);
	%Preprocessing(&mu_vars,'b0'  ,xb,b);
	%Preprocessing(&phi_vars,'d0'  ,wd,d);
	%Preprocessing(&zero_vars,'zero0',zeroxb,zero);
	%Preprocessing(&one_vars,'one0',onexb,one);

	proc nlmixed data=&Dataset tech=&tech &details;
		pizero=exp(&zeroxb)/(1+exp(&zeroxb));
		pione=exp(&onexb)/(1+exp(&onexb));
		mu=exp(&xb)/(1+exp(&xb));	
		phi=exp(&wd);
		w=mu*phi;
		t=phi - mu*phi;

		if(&depvar=0) then ll=log(pizero);
			else if(&depvar=1) then ll=log(pione);
			else ll=lgamma(w+t) - lgamma(w) - lgamma(t) + ((w-1)*log(&depvar)) + 
				((t-1)*log(1 - &depvar)) + log(1-pizero) + log(1-pione);

		model &depvar ~ general(ll);
		predict mu out=mu_results(keep=&depvar pred);
		predict phi out=phi_results(keep=&depvar pred);
		predict pizero out=pizero_results(keep=&depvar pred);
		predict pione out=pione_results(keep=&depvar pred);
	run;
	%Postprocessing(mu_hat,mu_results,&depvar);
	%Postprocessing(phi_hat,phi_results,&depvar);
	%Postprocessing(pizero_hat,pizero_results,&depvar);
	%Postprocessing(pione_hat,pione_results,&depvar);

	data prediction;
		merge mu_results phi_results pizero_results pione_results;
		by record;
	run;
%mend;












*Beta_Regression(
	Dataset,
	tech,
	details,
	mu_vars,
	phi_vars,
	zero_vars,
	one_vars,depvar);




%Beta_Regression(
	mydata,
	trureg,
	,
	stimulus timepoint ,
	stimulus timepoint,
	stimulus timepoint,
	stimulus timepoint,
	pct_sniff);


%Beta_Regression(
	mydata,
	trureg,
	,
	stimulus timepoint ,
	stimulus timepoint,
	stimulus timepoint,
	stimulus timepoint,
	pct_sniff);

%Beta_Regression(
	mydata,
	trureg,
	,
	food lion ,
	food lion,
	food lion,
	food lion,
	pct_sniff);

proc freq data=mydata;
table age breed; run;

proc genmod data=mydata;
class treatment (ref='cow') id timepoint (ref='2') breed;
model sniff_dur/visible = treatment   |sex;
repeated subject=id;
run;

proc genmod data=mydata;
class treatment (ref='cow') id timepoint (ref='2') breed;
model Ears_forw/visible = treatment  age_class ;
repeated subject=id;
run;

proc genmod data=mydata;
class treatment (ref='cow') id timepoint (ref='2') breed;
model sniffpc = treatment   ;
repeated subject=id;
run;


proc genmod data=mydata;
class treatment (ref='cow') id timepoint (ref='2') breed;
model approach_lat/visible = treatment  | age_class ;
repeated subject=id;
run;

/**/

proc glimmix data=mydata;
class treatment (ref='cow') id timepoint (ref='2');
model sniff_dur/visible = treatment | age_class treatment*sex
	/ solution ;
random _residual_ / subject=id type=cs; 
run;
proc glimmix data=mydata ;
class treatment (ref='cow') id timepoint (ref='2');
model sniff_dur/visible = treatment
	/ solution ;
random intercept / subject=id type=cs; 
run;

proc glimmix data=mydata;
class treatment (ref='cow') id timepoint (ref='2');
model Ears_forw/visible = treatment | age_class treatment*sex
	/ solution ;
random _residual_ / subject=id type=vc; 
run;

proc glimmix data=mydata;
class treatment (ref='cow') id timepoint (ref='2');
model Withdr_dur/visible = treatment | age_class treatment*sex
	/ solution ;
random  _residual_ / subject=id type=un; 
run;


/*====*/
proc genmod data=mydata;
class treatment (ref='cow') id timepoint (ref='2') breed;
model Ears_forw/visible = treatment | age_class treatment*sex  ;
repeated subject=id;
run;

proc genmod data=mydata;
class treatment (ref='cow') id timepoint (ref='2') breed;
model Withdr_dur/visible = treatment | age_class treatment*sex;
repeated subject=id;
run;


proc genmod data=mydata;
class treatment (ref='cow') id timepoint (ref='2') breed;
model Chew_dur/visible = treatment | age_class treatment*sex;
repeated subject=id;
run;

proc genmod data=mydata;
class treatment (ref='cow') id timepoint (ref='2') breed;
model sniff_dur/visible = treatment | age_class treatment*sex;
repeated subject=id;
run;


proc genmod data=mydata;
class treatment (ref='cow') id timepoint (ref='2') breed;
model Interact_dur/visible = treatment | age_class treatment*sex   ;
repeated subject=id;
run;


proc genmod data=mydata;
class treatment (ref='cow') id timepoint (ref='2') breed;
model latency_corrected/visible = treatment    ;
repeated subject=id;
run;











proc logistic data=mydata;
class  timepoint (ref='2') treatment (ref='food');
model sniff_dur/visible =  timepoint | treatment /expb ;
run;

proc logistic data=mydata;
class stimulus (ref='1') timepoint (ref='2');
model chewpc/visible = stimulus | timepoint /expb ;
run;


proc logistic data=mydata;
class stimulus (ref='1') timepoint (ref='2');
model earspc/visible = stimulus | timepoint /expb ;
run;



proc logistic data=mydata;
class stimulus (ref='1') timepoint (ref='2');
model withdrpc/visible = stimulus | timepoint /expb ;
run;


proc logistic data=mydata;
class stimulus (ref='1') timepoint (ref='2');
model interactpc/visible = stimulus  | timepoint /expb ;
run;
 
proc univariate data=mydata;
histogram;
var pct_sniff pct_chew pct_ears pct_withdraw pct_interact blows_occur;
run;

* This is the program for joint Analysis of Correlated Repeated Measures and Recurrent Events Processes 
* in the Presence of Death;

data one;	/* Repeated measures data */
input id y trt;
aa=2;
cards;
1 0.5 1 
1 1.2 1 
1 2.0 1 
1 3.4 1 
2 2.3 0 
2 2.5 0 
2 2.7 0 
...
;
run;
/*
 Status=1: Recurrent event; Event=2: death; Event=0: indepdent censoring;
*/
data two;
input id stoptime event trt;
aa=1;
cards;
1 5.0 1 1	/* recurrent event at time 5.0 */
1 6.7 1 1	/* recurrent event at time 6.7 */
1 10.0 2 1	/* death event at time 10.0 */
2 3.0 1 0
2 5.4 1 0
2 7.9 0 0	/* independent censoring event at time 7.9 */
3 4.5 1 1
3 6.0 0 1
......
;
run;
/* Data for Recurrent Events */
data recur;
set two;
if event=1;
run;
* Get the quantile for recurrent events;
proc univariate data=recur noprint;
var stoptime; 
output out=quant_r pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qr; 
run;
data quant_r;
set quant_r;
aa=1;
run;
/* Data for Death */
data death;
set two;
if event=2;
run;
* Get the quantile for death;
proc univariate data=death noprint;
var stoptime; 
output out=quant_d pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=qd; 
run;
data quant_d;
set quant_d;
aa=1;
run;
* Merge data with the quantiles;
data three;
merge two quant_r quant_d;
by aa;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data four;
set three;
array quant_r {11} qr0 qr10 qr20 qr30 qr40 qr50 qr60 qr70 qr80 qr90 qr100;
array quant_d {11} qd0 qd10 qd20 qd30 qd40 qd50 qd60 qd70 qd80 qd90 qd100;

array dur_r {10} dur_r1-dur_r10;
array dur_d {10} dur_d1-dur_d10;

array event_r {10} event_r1-event_r10;
array event_d {10} event_d1-event_d10;

do i=1 to 10;
	dur_r{i}=0;
	dur_d{i}=0;
end;

do i=1 to 10;
	event_r{i}=0;
	event_d{i}=0;
end;

* For recurrent event;
if event=1 then do;
	do i=2 to 11;
		if stoptime<=quant_r{i} then do;
			event_r{i-1}=1;		/* indicator of recurrent event in each interval */
			i=11;
		end;
	end;
end;

else do; /* If death or censored observation */
	do i=2 to 11;
		if stoptime<=quant_r{i} then do;	/* duration in each recurrent event quantiles */
			dur_r{i-1}=max(stoptime-quant_r{i-1}, 0);
			i=11;
		end;
		else dur_r{i-1}=quant_r{i}-quant_r{i-1};
	end;

	do i=2 to 11;	/* duration in each death event quantiles */
		if stoptime<=quant_d{i} then do;
			event_d{i-1}=(event=2);
			dur_d{i-1}=max(0, stoptime-quant_d{i-1}); 	/* indicator of death event in each interval */
			i=11;
		end;
		else dur_d{i-1}=quant_d{i}-quant_d{i-1};
	end;
end;

run;
data five;
set four one;
run;
data six;
set five;
if stoptime=. then stoptime=y;
run;
proc nlmixed data=six qpoints=5;

parms r01=0.1 r02=0.2 r03=0.3 r04=0.4 r05=0.5 r06=0.6 r07=0.7 r08=0.8 r09=0.9 r10=1.0
		h01=0.025 h02=0.05 h03=0.075 h04=0.1 h05=0.125 	h06=0.15 h07=0.2 h08=0.25 h09=0.3 h10=0.4 
		beta1=1 eta1=1 alpha0=0 alpha1=1 gamma1=-1.5 gamma2=-.5 gamma3=1.5
		vara=1 varb=.5 vare=1;
bounds r01 r02 r03 r04 r05 r06 r07 r08 r09 r10 
h01 h02 h03 h04 h05 h06 h07 h08 h09 h10 vara varb >=0;

if aa=2 then do;	/* likelihood for repeated measures y */
	mu3= alpha0 +  alpha1 * trt +  a ;
	loglik=-.5*(y-mu3)**2/vare-.5*log(2*3.14159*vare);
end;

if aa=1 then do;	/* likelihood for recurrent events and death */

* base hazard and cumulative baseline hazard for recurrent events;

base_haz_r=r01 * event_r1 + r02 * event_r2 + r03 * event_r3 + r04 * event_r4 + r05 * event_r5 + 
r06 * event_r6 + r07 * event_r7 + r08* event_r8 +r09 * event_r9 + r10 * event_r10;
cum_base_haz_r=r01 * dur_r1 + r02 * dur_r2 + r03 * dur_r3 + r04 * dur_r4 + r05 * dur_r5 + r06 * dur_r6 + 
r07 * dur_r7 + r08* dur_r8 +r09 * dur_r9 + r10 * dur_r10;

* base hazard and cumulative baseline hazard for death events;

base_haz_d=h01 * event_d1 + h02 * event_d2 + h03 * event_d3 + h04 * event_d4 + h05 * event_d5 + 
h06 * event_d6 + h07 * event_d7 + h08* event_d8 +h09 * event_d9 + h10 * event_d10;
cum_base_haz_d=h01 * dur_d1 + h02 * dur_d2 + h03 * dur_d3 + h04 * dur_d4 + h05 * dur_d5 + h06 * dur_d6 + 
h07 * dur_d7 + h08* dur_d8 +h09 * dur_d9 + h10 * dur_d10;

mu1= beta1 * trt + gamma1 * a + b;			/* for recurrent event */

mu2= eta1 * trt  + gamma2 * a + gamma3 * b;	/* for death event */

loglik1=-exp(mu1) * cum_base_haz_r;

loglik2=-exp(mu2) * cum_base_haz_d;

if event=1 then loglik= log(base_haz_r) + mu1 ; 					/*log likelihood for recurrent event */
if event=2 then loglik=loglik1 + log(base_haz_d) + mu2 + loglik2;	/*log likelihood for death */
if event=0 then loglik=loglik1 + loglik2;							/*log likelihood for censoring */

end; /* End if aa=1 then do */	

model stoptime ~ general(loglik);
random a b ~ normal([0, 0], [vara, 0, varb]) subject=id;
ods output ParameterEstimates=est1_all FitStatistics=fit1_all; 
run;

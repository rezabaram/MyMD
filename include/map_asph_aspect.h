#ifndef MAP_ASPH_ASPECT_H
#define MAP_ASPH_ASPECT_H 
#include"exception.h"

//mapping between asphericity and aspect ratio for oblate and prolate spheroids 
//brute force :-)

double asph_start=1.005;
double asph_step=0.005;
double asph_end=1.2;
double asph[40]={1.005,
	1.01,
	1.015,
	1.02,
	1.025,
	1.03,
	1.035,
	1.04,
	1.045,
	1.05,
	1.055,
	1.06,
	1.065,
	1.07,
	1.075,
	1.08,
	1.085,
	1.09,
	1.095,
	1.1,
	1.105,
	1.11,
	1.115,
	1.12,
	1.125,
	1.13,
	1.135,
	1.14,
	1.145,
	1.15,
	1.155,
	1.16,
	1.165,
	1.17,
	1.175,
	1.18,
	1.185,
	1.19,
	1.195,
	1.2};


double aspect_oblate[40]={0.8474612999057303,
0.7923027162407834,
0.7527887576668176,
0.7212430308590593,
0.6947141820788358,
0.6716991492128184,
0.6513110901228971,
0.6329760067141943,
0.6162977948850417,
0.6009899038611373,
0.5868373802344924,
0.5736742706336806,
0.5613694166079694,
0.5498171253543555,
0.5389308176680306,
0.5286385701792486,
0.5188799053641128,
0.5096034281835522,
0.5007650520788103,
0.49232664452099373,
0.48425497718940896,
0.47652090124740537,
0.469098691578197,
0.46196551964662097,
0.45510102554030485,
0.4484869673784151,
0.44210693171673876,
0.4359460925121367,
0.4299910090928861,
0.4242294557210768,
0.4186502769390005,
0.4132432641095829,
0.40799904949370264,
0.40290901492840775,
0.3979652127321625,
0.39316029690494186,
0.38848746304069615,
0.38394039564845717,
0.37951322180208036,
0.37520047021924724};

double aspect_prolate[40]={1.185273387204022,
1.2734586036783162,
1.3462928551363025,
1.4114542492582445,
1.4718958540749887,
1.5291224458861417,
1.584022580868675,
1.6371719709685262,
1.688968433911266,
1.7397002338297263,
1.7895840371878236,
1.8387875100741922,
1.8874435232829732,
1.9356594815526957,
1.9835236755827808,
2.0311097397646467,
2.0784798621398877,
2.1256871477312513,
2.1727773925180704,
2.219790437858601,
2.2667612202845326,
2.313720596198884,
2.3606959976148856,
2.4077119592698177,
2.454790546560319,
2.5019517061110004,
2.5492135553475954,
2.5965926235116,
2.6441040536698135,
2.69176177313272,
2.7395786380896103,
2.787566557050489,
2.835736596751893,
2.8840990734625334,
2.9326636320626576,
2.981439314829362,
3.030434621510161,
3.0796575619886837,
3.1291157026223995,
3.1788162071517596};

double get_aspect_oblate( double a){
	ERROR(a<asph_start or a>asph_end,"Out of bound asphericity requested. The valid range is: [1.005,1.2]");
	int i=(a-asph_start+10e-10)/asph_step;
	return aspect_oblate[i];
	}

double get_aspect_prolate( double a){
	ERROR(a<asph_start or a>asph_end,"Out of bound asphericity requested. The valid range is: [1.005,1.2]");
	int i=(a-asph_start+10e-10)/asph_step;
	return aspect_prolate[i];
	}



#endif /* MAP_ASPH_ASPECT_H */

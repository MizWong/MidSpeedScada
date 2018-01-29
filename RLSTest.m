clear all;

clear fnRLS;

x  = [-0.886120723168321,0.192680685592617,0.296878614985397,0.0261986172176128,-0.0297580662552420,-1.79978033395490,-0.284966266391269,1.47564107985029,-0.695115941520914,0.447216644379263,0.886604865018915,1.11405915005748,2.01939101789571,0.0107139669724275,-0.981906337883993,-0.275275923250526,-1.18152689144476,0.361641274297749,0.920345260176028,-3.07448182247052,0.526065566115557,1.42683463203719,0.970956283747279,1.41884841867172,0.618367198714028,-0.492348580237500,-0.347920124149852,0.0703524739784332,0.0293008776413065,-1.48135521523594,0.822807436115081,-1.66437897063440,-0.654024241123598,-2.16730098780970,0.585829762477305,0.0457877274662108,-0.845638392246517,0.116931603356211,-1.02240110869994,-0.618099496813119,0.475186989117375,1.09082471618658,-0.106338154172232,1.09501311174362,0.834865037997992,1.79782108351462,-1.29450055276691,0.451565106807863,0.768263115221266,-0.526535539193169,0.900661240528551,0.690109140286225,0.0451273044309528,-1.13025874092724,-1.36655244368176,0.378250470812870,-2.23903910335473,1.67720786904321,-0.251942341960518,-0.319515109433189,-0.858307803718376,-0.695739500549735,-0.494581125670444,-0.0996958453584543,0.853537697516668,1.06590878356558,1.89972229162429,-0.404002259724751,-0.691249111045340,0.0294409235048039,0.473468705633933,0.995885916187622,-0.614941565505588,0.773476898279406,-0.682193848693788,-0.486417819138514,-1.75051293790162,-0.0857001796291393,0.704778187695878,-0.183651311780593,-0.278882588224000,0.127947352519338,-0.959286267736922,0.158256187012633,-1.13969192347908,-0.681869105441182,0.0711016607532895,0.876121687557969,0.661135748305552,1.17680087169485,0.699757115669452,2.01900645748326,-0.385506740721730,-1.23352705176606,-0.0725490936754995,0.0220997269582148,0.378819504311197,0.199136410319910,-1.77757645294406,-2.05268829698405,-0.0267700472261589,-0.224864576366999,-0.224853081038242,0.128205722833401,0.545311505022903,1.66074383541494,0.741565604202386,0.0573541161769726,0.127822363800496,1.08146883290683,0.962265318599557,0.794514528594131,-1.48017651161543,0.874897392459751,0.719544681677914,0.818557270693202,-1.66417983422068,0.00632907180909204,0.0740481109423326,-0.881326902366911,0.668932403770361,-0.136603898147037,-0.519336312240710,0.900279936817984,-1.64846858448884,0.724818015409220,-0.501079522468061,0.643625496958633,-0.241948633843596,0.797720308275797,-0.277767891528426,-1.68119143893146,0.484054493997490,-0.0355662900278660,0.130006484917379,-0.395602277652553,-1.08366342984429,-0.651315880369764,0.317893600843155,1.25122002329851,-0.746668167728688,-0.462300680842573,-1.58473090345680,1.71038216407685,-0.804987558018463,0.818210404335906,-0.0262590925167039,-2.10604216323917,0.350787855903339,0.962432539063534,0.763406945645507,0.798557360475981,-2.31139454767724,0.373843224492469,1.26986426363615,0.622711344923395,-1.76509038556298,0.801732647494845,-0.550031940412051,-0.460500501072916,0.319798534415112,-0.852450651283526,0.0419003637349001,0.0610688722210733,-0.800091335109284,0.574684133223634,1.31202821373708,-0.149247363318893,1.48128402206159,-0.281279027755726,-1.13578679851025,-0.902132582140538,1.74589461394780,1.09294823804322,-0.243781167397389,1.79314067582924,-0.404668527425120,2.62506117290121,0.174827910115855,2.24702660150789,-0.243710396088646,0.625254523830055,-0.436429585297289,1.54345343812076,-1.05837065520087,0.397745062891466,-1.34065941411646,0.0321183573239437,-0.505680082632591,0.542293908718360,-0.407240861314871,0.175113182133515,0.249472925087296,-0.322863395566212,0.734465322999983,0.518610349287192,-1.13204665153937,-2.05191765125932,-0.367349533110419,0.391573533495298,-0.269966102087105,-1.94337200343129,0.549977389606579,-0.666132336242961,-1.05157501646622,0.888419992243212,-0.531580581247482,0.455898866201915,0.235100534808325,1.27023868082095,1.61768281556414,-0.0713826439855842,0.823430076675969,0.692336586007272,-0.385697610619369,-0.428682544553800,1.22085103793882,1.36758642052428,1.93679381562407,0.0978109652020670,0.438088572307358,1.01559136929510,-1.05087106525105,-0.171531420460394,0.498409606547821,1.11448121639513,0.329983070022336,0.231602796610206,1.08532281403169,-1.04236582767265,-0.357517481821262,1.97360834221646,-0.168880355515087,0.553014769074983,0.698502966241002,-0.228559544772605,-1.18213427705137,-0.905416533736302,0.755867455533863,-1.11577035412283,-1.95875258799419,1.24159786345581,0.0879963564257252,-0.462585654140075,-0.263218665585316,-0.126145324865788,0.0997727802266499,0.522666332924460,-0.351562394296747,0.531386147094366,-0.281329292500329,1.19027647570129,-0.0773487346129768,1.43390476202202,-1.38109333690480,0.299587065535568,-0.456547837263391,0.415282653124032,0.691606189028949,-0.306453911086789,-0.719092565245578,0.601925612295982,0.424751866462079,-0.0103916733892786,-0.588389945316967,-0.960333511548949,1.55962291820452,-0.707325932559535,-0.286084990081098,1.25517444354342,1.30367554419848,-1.53307105820794,0.392645519698888,-0.447197546419474,0.144304091457626,-0.893234426670985,0.478477735424532,-0.149079945803450,-1.97932534419626,-1.47243202390925,0.860415634505840,0.658720558946802,1.55373594446208,-0.671499417925613,0.508716009921796,1.54212703817771,-2.01235791880081,0.731494254334072,0.323854960289336,-1.18493012408949,0.0636224584837511,-1.09117750286005,-1.01575839782779,-0.726717987336530,1.03269390606190,-0.323643502121636,0.747928731170256,-1.68394240918828,0.309169981957823,-0.654694668233158,0.215059665630235,-0.371924881404312,-0.262975218526364,-0.144133872579392,-1.71005936316219,1.26554047168748,0.315500989899835,-0.0208468823114189,0.243063122284456,1.33731092006264,1.60218789542812,0.614443981954428,-0.0335659432364702,-0.198284122101138,-1.13579782125012,-1.34149149541688,-1.26826093906505,0.411999365374415,-0.641128573076002,-0.815621783662883,0.817574763707832,1.54243442864750,-0.322326770238622,-0.185693394792511,-1.93064321671421,2.07449869911780,-0.358604951683831,-0.913034212151687,-1.79197132280219,-0.408085549968542,0.719362768343945,0.887602049526791,0.238564494163894,0.665475142629064,-1.00041828921855,1.14509937285861,-0.287802345339847,-0.250593339302539,-0.886921657795571,-0.479886202564105,-0.957029518674574,-0.0162927641498263,0.891226353595894,1.36938503602519,2.61971216641514,-0.584499269730205,-0.849191457596360,0.0631926391363570,0.484290258051757,1.85224076686872,0.0671253777435237,-1.90391984439207,0.388680159915724,-1.34389236616605,-0.186604965300434,0.998838691922671,0.732074550198000,-0.0162690795113025,0.707132387344537,-0.479558095559357,-1.27243567053652,-0.179344460743519,0.200251820619224,0.189706867597281,-0.372721041280032,-1.12342556589379,-1.95812701221574,-1.15166218810205,-0.112811396335503,0.189545749361032,1.41634964323008,-0.568282605270808,-0.645214292190800,0.341164550884576,-0.677452843607989,0.100990396886530,-0.919716686808636,0.328622826665100,0.519304376031693,-0.274030014490412,0.0907915356559087,1.05909356936814,1.47670863516924,-0.247999418034607,-0.522702097863536,-1.58691732457943,-0.405781661362432,0.326426438070460,0.0709300453613238,1.00862382038353,0.596816326513099,-0.884938787951747,-0.537975527012474,-0.644254739421211,-1.11849769101592,-0.619847382506099,0.819095343763675,-0.734939209803989,0.319958975427265,-1.28367690047790,1.75242746412388,0.0638249547912932,-2.97835671775048,-0.808862651773318,0.0564117507874823,-0.0471354103952045,-0.532651369297542,0.481786449203947,0.979438632866504,0.526293116070362,1.52067356074261,-1.65368638846988,-0.0581245021453148,1.37032757033428,1.24165481637056,-0.585708553691064,1.16322608241206,1.43606096355852,0.857925203916704,1.15435761685033,1.05941251465803,-0.0306361234412506,-0.0785481164793454,1.45700253636412,0.280362870619959,-0.343445680709679,0.123104009494519,-1.09309378938023,0.402432473090381,-0.948431087024061,-1.25038269539271,-0.201098190260943,-1.25255135232683,2.60391416422096,-1.15575471539870,0.914242062082552,0.434540299712368,1.26035197170137,-2.75527906594141,-2.24536766025217,0.667902382716336,-0.886534253179635,-0.630196816688762,0.344513672879003,0.762197715930989,1.22588980228870,2.06331886522023,0.821437095049616,-0.704638790848611,1.09769040252077,-0.367771982586970,1.24925858623572,-0.0488586237277714,-0.466725274001515,0.304249217999359,-1.64323107365409,0.0822533214884235,0.0700480712618135,-1.13835960796258,0.00863537635127596,-0.529457706114655,0.0202464800684006,-0.565690553425667,0.121143519236406,1.68172613548000,0.0901968488726735,-0.327998315296923,0.251287598114421,0.706627158573130,-0.201582942506191,0.467168387201280,-1.40557868726787,-0.734867901292414,0.358467781896503,0.0242329487065203,-1.02080521841007,0.517686679036942,-1.17596754637588,-0.188102951638436,0.0375765673791731,1.00437656988688,0.432840368097017,1.35806094971474,0.174221712856769,-1.24757255605820,0.0224384288437203,1.06618408675834,-0.146863351865547,-1.46419491270580,-0.654510258873041,-1.51278994949867,1.32636033974428,1.06843885941372,-0.275605940014035,-1.25321474077060,1.07301163948530,-0.290196385985675,-1.35043246362408,0.347584443541381,-0.533214720918366];
b  = fir1(31,0.5);     % FIR system to be identified
n  = 0; % Observation noise signal
d  = filter(b,1,x)+n;  % Desired signal
P0 = 10*eye(32); % Initial sqrt correlation matrix inverse
lam = 0.99;            % RLS forgetting factor
% ha = adaptfilt.rls(32,lam,P0);
% [y,e] = filter(ha,x,d);
% plot(1:500,[y]);
% % xlabel('Time Index'); ylabel('Signal Value');
% % subplot(2,1,2); stem([b.',ha.Coefficients.']);
% % legend('Actual','Estimated');
% % xlabel('Coefficient #'); ylabel('Coefficient valUe');  grid on;

w = zeros(32,1);
y = zeros(1,length(x));

hRLS = fnRLSCreate(32, lam);

for ii = 1:1:length(x)
    if ii < 33
        y(ii) = [zeros(1,32-ii),x(1:ii)] * w;
        e = d(ii) - y(ii);
        [deltaW, hRLS] = fnRLS( hRLS, [zeros(1,32-ii),x(1:ii)]', e, 0 );
        w = w + deltaW;
    else
        y(ii) = x(ii-32+1:ii) * w;
        e = d(ii) - y(ii);
        [deltaW, hRLS] = fnRLS( hRLS, x(ii-32+1:ii)', e, 0 );
        w = w + deltaW;
    end
end


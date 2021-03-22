(* ::Package:: *)

(* ::Title:: *)
(*Deflector Data Processing*)


(* ::Chapter:: *)
(*2021-03-22*)





BeginPackage["deflectorProcessing`"];

PlotDeflectorData::usage="PlotDeflectorData[data,scale,plotLabel] Plots data of the form {{v1_1,v1_2,sig1},{v2_1,v2_2,sig2},...}, in an array plot";
OrganizeDeflectorData::usage="OrganizeDeflectorData[filename,voltagePairs] Returns voltage pairs list with third item, the current";
DataExclude::usage="DataExclude[dat_,eastLowerBound_,westLowerBound_]";


Begin["`Private`"];
DataExclude[dat_,eastLowerBound_,westLowerBound_]:=Select[dat,Abs[#[[1]]]>=eastLowerBound && Abs[#[[2]]] >= westLowerBound&]

OrganizeDeflectorData[fn_(*filename*),vp_(*voltage pairs*)]:=Module[{sr=11(*startRow*),
er=-2(*endRow*),
ri=Import[fn](* raw import*),
v (*values*), 
sig(*signal*),
graphs},
v=ri[[sr;;er]]; (* (v)alues *)
sig=Transpose[v];(* There are columns for error, but currently they aren't populated. So we only take the odd columns *)
sig=<|"fd"->sig[[1]],(*fd=faraday collector*)
"ct"->sig[[3]],(* ct=chiral top *)
"cb"->sig[[5]], (* cb=chiral bottom *)
"he"->sig[[7]]|>; (* he=helium target*)
sig["he"]=sig["he"]*-1*^-9;


sig=Map[Flatten[Riffle[vp,#]]&,sig]; (* Combine signal and voltage pairs *)
sig= Map[Partition[#,3]&,sig] (* Organize the signal and voltage pairs into {east voltage, west voltage, signal} *)];

Clear[PlotDeflectorData]
PlotDeflectorData[sig_(* signal *),sc_(* scale *),
pl_(*Plot Label*),minMax_:None,ft_:{True,True}(*frame ticks*)]:=Module[
{vEast,vWest,temp,tempSig,
vMax(*voltage max*),
vMin(*voltage min *),
arr(* The list of arrays of values we ultimately plot *),
m (* a single array (matrix) *),s (* The signal in the for loop*),
i,j(* iteration variables *),
mm (* The min and max of the signal *),
me (* The mantissa and exponent of the maximum value.*),
scale,(*factor to scale the data by so the legend can bemore compact.*)
colorFunc,viridisColors,
graphs},

(* The quirk of the array plot is that its for integers, so we have to convert our array of values to an integer so that they can fit into a table. *)
(* It's much easier to think about rescaling data when the max is the furthest number away from 0 and the min is the closest number to zero, so I make our negative currents into positive. *)
tempSig=Transpose[{Transpose[sig][[1]],Transpose[sig][[2]],-Transpose[sig][[3]]}];
temp=Transpose[tempSig];
vMax=Max[Catenate[temp[[;;2]]]];
vMin=Min[Catenate[temp[[;;2]]]];
mm=MinMax[temp[[3]]];
If[mm[[1]]<0,mm[[1]]=1*^-16];
(* These lines are for the better bar legend which has the exponent above the bar *)
me=MantissaExponent[mm[[2]]]*{10,1}+{0,-1}; (* Move the mantissa to be a value between 10 and 1/10, shift the exponent accordingly *)
scale=Power[10,me[[2]]];
mm=mm/scale;

viridisColors={RGBColor[0.2823645529290169, 0., 0.3310101940118055],RGBColor[0.28442016331743475`, 0., 0.3370344273837448],RGBColor[0.28640872215228885`, 0., 0.34299986027156143`],RGBColor[0.28832801827269444`, 0., 0.3489014246790825],RGBColor[0.29017937839961166`, 0.0026548744457282002`, 0.3547381197779056],RGBColor[0.2919602785108896, 0.007927574441840031, 0.3605089628592038],RGBColor[0.2936725695568318, 0.013523701340444767`, 0.3662108687467896],RGBColor[0.2953145929546483, 0.0194469283736754, 0.37184193824445494`],RGBColor[0.2968856990777091, 0.025709272843022362`, 0.377402604549329],RGBColor[0.2983864384776219, 0.03232004345340847, 0.38288925650020894`],RGBColor[0.29981597062377513`, 0.039284989212368455`, 0.38829907434247934`],RGBColor[0.3011730649474346, 0.04626360509301145, 0.39363381249817775`],RGBColor[0.3024574616405568, 0.05292448251262403, 0.39888773777233144`],RGBColor[0.3036695909467026, 0.059347679387570774`, 0.4040625134864877],RGBColor[0.3048077646670268, 0.0655767750733681, 0.4091545009758388],RGBColor[0.3058726942921675, 0.07164698633434458, 0.4141620022087494],RGBColor[0.30686354185634795`, 0.07758748304703458, 0.419085269399855],RGBColor[0.3077802784412034, 0.0834178665180236, 0.4239197611787469],RGBColor[0.3086224941211352, 0.08915477951604339, 0.4286679531263877],RGBColor[0.3093909368237089, 0.09481086300337191, 0.43332440286164786`],RGBColor[0.3100840531043614, 0.10039934660784217`, 0.4378904725411851],RGBColor[0.3107021547526659, 0.10592544707225837`, 0.4423628376506268],RGBColor[0.3112451783861907, 0.111398995007257, 0.44674164765622],RGBColor[0.3117137300318166, 0.11682684097453058`, 0.45102599429604956`],RGBColor[0.31210749532005316`, 0.12221105850499253`, 0.45521366881321973`],RGBColor[0.3124264208600718, 0.12756000041340385`, 0.45930580962347606`],RGBColor[0.31267102574316064`, 0.1328717148189205, 0.4632993229731653],RGBColor[0.31284284828449566`, 0.13815390993203355`, 0.46719523432267124`],RGBColor[0.31293985716270956`, 0.1434064187456714, 0.4709905542704291],RGBColor[0.3129647555680025, 0.14863302233253345`, 0.4746882091147527],RGBColor[0.3129172762828141, 0.15383324394186293`, 0.47828532474775726`],RGBColor[0.3127980600772598, 0.1590082922293873, 0.48178183314421624`],RGBColor[0.3126079674503442, 0.16416130934378315`, 0.4851796030849527],RGBColor[0.3123476422504348, 0.16929153088387608`, 0.48847588346394805`],RGBColor[0.3120187867754984, 0.17439996402117275`, 0.4916733669375038],RGBColor[0.3116232941035938, 0.17948560856731854`, 0.4947714272765002],RGBColor[0.31115984744015995`, 0.1845516383459222, 0.4977706233250952],RGBColor[0.3106312638671459, 0.189595868940358, 0.5006704469605663],RGBColor[0.31003859809087064`, 0.19461776801237182`, 0.5034736673089839],RGBColor[0.30938382493792904`, 0.19962012298531123`, 0.5061799701305794],RGBColor[0.3086683938758612, 0.20459957100719467`, 0.5087899731598793],RGBColor[0.30789302647936884`, 0.20955888696372615`, 0.5113070173032073],RGBColor[0.307058669956948, 0.21449483201644065`, 0.513728818335263],RGBColor[0.30616963655061763`, 0.2194100154940718, 0.5160606205820103],RGBColor[0.3052269065043157, 0.22430118261968157`, 0.5183002465878732],RGBColor[0.3042309077920497, 0.22916966474571626`, 0.520452297156245],RGBColor[0.30318432730096756`, 0.2340159100780761, 0.5225163029836838],RGBColor[0.30208910989145576`, 0.23883855539609794`, 0.5244952101282286],RGBColor[0.3009472750064046, 0.24363810897601756`, 0.5263910452541571],RGBColor[0.29976132132321154`, 0.24841118674175394`, 0.5282048535123068],RGBColor[0.29853340919581994`, 0.2531623990240717, 0.5299405582099276],RGBColor[0.29726490390425825`, 0.2578873976329164, 0.5315972918179168],RGBColor[0.29595755856435807`, 0.2625874493056926, 0.5331785030475249],RGBColor[0.2946138805993375, 0.26726378849777516`, 0.5346875949663166],RGBColor[0.2932353894530138, 0.27191415168129185`, 0.5361247892898974],RGBColor[0.2918258704725537, 0.27653989400666285`, 0.5374938454348637],RGBColor[0.290387139948913, 0.28113768475292045`, 0.5387960766501153],RGBColor[0.2889195372181988, 0.28571303833310896`, 0.5400351419012158],RGBColor[0.28742616289551665`, 0.29026061879591203`, 0.541210431546397],RGBColor[0.28590907506551183`, 0.29478258938071067`, 0.5423272767396595],RGBColor[0.28436845691665213`, 0.2992802217180254, 0.5433859495199381],RGBColor[0.28280885381705473`, 0.30375224020966846`, 0.5443899159864241],RGBColor[0.28123118105327793`, 0.30819892495661744`, 0.5453407090892604],RGBColor[0.2796382036694848, 0.3126190607542608, 0.5462408516171586],RGBColor[0.2780293659165648, 0.3170149456432699, 0.5470928102383564],RGBColor[0.27640874852219066`, 0.321385414500347, 0.5478981611011682],RGBColor[0.27477620046950885`, 0.3257316254417039, 0.5486591285070745],RGBColor[0.2731336070036193, 0.3300547425253127, 0.5493769038436812],RGBColor[0.2714830730483761, 0.3343525579890229, 0.5500541290065128],RGBColor[0.26982545321488427`, 0.33862832893130074`, 0.550693159876779],RGBColor[0.2681626147651372, 0.3428799001623282, 0.5512956789400659],RGBColor[0.26649507924006577`, 0.34710848826302004`, 0.5518619993206375],RGBColor[0.26482494588026595`, 0.35131500140549954`, 0.5523958384505545],RGBColor[0.2631526559675265, 0.3554995323434884, 0.5528983516354158],RGBColor[0.261479290708243, 0.3596631879129328, 0.553370667054746],RGBColor[0.2598068996301584, 0.36380486101494103`, 0.5538155494173054],RGBColor[0.2581351267096509, 0.36792673591505964`, 0.5542332292091459],RGBColor[0.2564654875671387, 0.3720277792587921, 0.5546244907783647],RGBColor[0.25479666350991237`, 0.3761100924489166, 0.5549915450448365],RGBColor[0.25313154461533993`, 0.3801736871524052, 0.5553351998790793],RGBColor[0.25146893469049797`, 0.384219609636479, 0.5556576460945991],RGBColor[0.24980931862578432`, 0.3882468863669942, 0.5559587153582112],RGBColor[0.24815451941163053`, 0.39225657410893083`, 0.5562394960910725],RGBColor[0.24650281919380035`, 0.3962507215819094, 0.5565020591903272],RGBColor[0.2448554379333428, 0.40022836834946646`, 0.556746262161579],RGBColor[0.24321127505789647`, 0.40419159051679776`, 0.5569722290778042],RGBColor[0.24157208508670241`, 0.40813938952902734`, 0.557181850680197],RGBColor[0.23993532100520137`, 0.4120738222862286, 0.5573752305636125],RGBColor[0.23830244481673896`, 0.41599392773534455`, 0.5575532711989851],RGBColor[0.23667266390977104`, 0.41990177352141406`, 0.5577160128723664],RGBColor[0.23504468776110826`, 0.4237973655416275, 0.5578634992121222],RGBColor[0.23341840155293508`, 0.42768077962221346`, 0.5579966435791459],RGBColor[0.23179414987893965`, 0.43155402955947875`, 0.5581155142003147],RGBColor[0.23017061991847487`, 0.4354151520581061, 0.5582210548062131],RGBColor[0.2285458366653788, 0.4392682409766716, 0.5583112863650675],RGBColor[0.2269217957032748, 0.44311029176723954`, 0.5583881853088769],RGBColor[0.2252950588660499, 0.4469443851124765, 0.5584507478170652],RGBColor[0.22366535734333398`, 0.45077047699087197`, 0.5584999999185336],RGBColor[0.22203317321082194`, 0.45458770680091215`, 0.5585338908615732],RGBColor[0.22039568222953204`, 0.4583990258215093, 0.5585544415094105],RGBColor[0.2187534764148667, 0.4622025551537203, 0.5585596247099263],RGBColor[0.21710379387686735`, 0.466001277474918, 0.5585494296782225],RGBColor[0.21544640014275412`, 0.46979327601523685`, 0.5585238597828509],RGBColor[0.2137812367677258, 0.4735795360665342, 0.5584829093200094],RGBColor[0.21210435581990783`, 0.4773631248343298, 0.5584245259010897],RGBColor[0.21041654636306228`, 0.4811400844455726, 0.5583487508089193],RGBColor[0.2087156715163057, 0.48491543980046187`, 0.5582545074731591],RGBColor[0.20700265458437167`, 0.4886862374146865, 0.5581428559806787],RGBColor[0.2052740577019755, 0.49245445289710926`, 0.5580117160916531],RGBColor[0.20352973151672893`, 0.4962192161230187, 0.5578601343029762],RGBColor[0.20176907507316538`, 0.49998151568117394`, 0.5576870630720852],RGBColor[0.19999051193019307`, 0.5037433207187862, 0.5574934462432263],RGBColor[0.1981942071395154, 0.50750175769624, 0.5572763581459699],RGBColor[0.19637893240486903`, 0.5112607523934385, 0.5570376811361709],RGBColor[0.19454390995867465`, 0.5150164307230379, 0.5567735144515659],RGBColor[0.1926877983020047, 0.5187727514082403, 0.5564836952783472],RGBColor[0.1908096988450276, 0.522527822045104, 0.556167356233582],RGBColor[0.1889101064832561, 0.5262826202118115, 0.5558234063897372],RGBColor[0.1869882250558732, 0.5300371054507527, 0.5554507440723478],RGBColor[0.18504404737541283`, 0.5337914192017983, 0.5550485162697486],RGBColor[0.18307650380221172`, 0.5375464935926217, 0.5546145021343457],RGBColor[0.18108751106392554`, 0.5413014311482989, 0.5541489011713485],RGBColor[0.179079069917763, 0.5450561095248925, 0.5536514862639709],RGBColor[0.17704851110868655`, 0.5488117187514125, 0.5531194414635993],RGBColor[0.17499810635812427`, 0.5525661768772355, 0.5525516502615303],RGBColor[0.1729271521674646, 0.5563235081702055, 0.551946905668582],RGBColor[0.17084034495202932`, 0.5600797579540071, 0.5513055188601756],RGBColor[0.16873518608371835`, 0.5638369045836091, 0.5506231120862802],RGBColor[0.16661789658301376`, 0.5675940102752802, 0.5499020260312123],RGBColor[0.16448693755461796`, 0.5713520340493227, 0.5491378466444521],RGBColor[0.16235010494692595`, 0.5751100369776962, 0.5483329604685923],RGBColor[0.1602069832547893, 0.5788679596972927, 0.5474829147739082],RGBColor[0.15806264119618174`, 0.5826258959527332, 0.5465891203401918],RGBColor[0.15592085037000353`, 0.5863828097619445, 0.5456483318142806],RGBColor[0.15378524571480254`, 0.5901406674305801, 0.5446602664371812],RGBColor[0.15166304622295707`, 0.593897573685252, 0.5436233831996301],RGBColor[0.14955934468428786`, 0.597654419423444, 0.5425371443011318],RGBColor[0.14748163623723387`, 0.6014093240653062, 0.5413990511580228],RGBColor[0.14543738134019324`, 0.6051641624269574, 0.5402105175168492],RGBColor[0.14343467030926343`, 0.6089160395607421, 0.5389691120794025],RGBColor[0.14148086657706913`, 0.6126659299137229, 0.5376724808099055],RGBColor[0.13958716457678366`, 0.6164157926020056, 0.5363202416485567],RGBColor[0.13776366739068555`, 0.6201616886553919, 0.5349120639299949],RGBColor[0.13602091941109576`, 0.6239065451618546, 0.5334461775790392],RGBColor[0.13437265977205576`, 0.6276464255764156, 0.5319223153351261],RGBColor[0.13283238318759788`, 0.6313842385203057, 0.5303396509309062],RGBColor[0.13141083773044995`, 0.6351181073443843, 0.528695928325246],RGBColor[0.13012562239437234`, 0.6388469415665682, 0.5269907342678339],RGBColor[0.12898908642655502`, 0.6425726931591681, 0.5252235696860159],RGBColor[0.12801762810243703`, 0.6462914463669078, 0.5233932929927063],RGBColor[0.1272255861040811, 0.6500071247171884, 0.5214978927023376],RGBColor[0.12663084317706338`, 0.6537147788902224, 0.519538336475508],RGBColor[0.12624680080960696`, 0.6574183247714932, 0.5175125345875667],RGBColor[0.12608994360406375`, 0.6611128336141637, 0.5154205123723917],RGBColor[0.12617263740036105`, 0.6648012725844235, 0.5132606594914964],RGBColor[0.1265110859590127, 0.6684815267484132, 0.511032405517675],RGBColor[0.12711466919306777`, 0.6721537530870352, 0.5087347666999795],RGBColor[0.12799581674473437`, 0.6758177719435436, 0.5063675657201724],RGBColor[0.12916500542967113`, 0.6794706989916038, 0.5039309536083489],RGBColor[0.13062688656598942`, 0.6831154175921678, 0.5014215766343267],RGBColor[0.13238792426640142`, 0.6867490315471613, 0.49884168959711156`],RGBColor[0.1344529141545566, 0.6903714630144, 0.4961885702691708],RGBColor[0.13682278693871636`, 0.6939836082650025, 0.493462471806926],RGBColor[0.13949652044856853`, 0.6975826171487558, 0.49066268626229187`],RGBColor[0.14247497810778012`, 0.7011702998442422, 0.48778873213754487`],RGBColor[0.14575022774163957`, 0.704742819719188, 0.4848389747066589],RGBColor[0.14932184918451274`, 0.7083029725794211, 0.48181384166285274`],RGBColor[0.1531793767768577, 0.7118479289119495, 0.47871379185797663`],RGBColor[0.15731856380201292`, 0.7153785887738957, 0.4755369438569324],RGBColor[0.16173071601073405`, 0.7188938119708388, 0.4722824057138442],RGBColor[0.16640600246038467`, 0.7223917868234861, 0.46895172451554185`],RGBColor[0.17133854153234118`, 0.7258742839343467, 0.4655431005241239],RGBColor[0.17651571973890343`, 0.7293375175127107, 0.46205513288486627`],RGBColor[0.1819300525085049, 0.7327852410807962, 0.45848792117067955`],RGBColor[0.1875710313130601, 0.7362116391655383, 0.45484323921432873`],RGBColor[0.19342964878541258`, 0.7396196144261882, 0.45111796445512237`],RGBColor[0.19949918559250934`, 0.7430089883586135, 0.44731304660604426`],RGBColor[0.2057678869700838, 0.7463760005250465, 0.4434292850786966],RGBColor[0.21222937521669094`, 0.7497233760041233, 0.4394635038573973],RGBColor[0.21887342785296, 0.7530473694983222, 0.4354175943073504],RGBColor[0.2256970102409276, 0.7563496686461588, 0.43128829183503475`],RGBColor[0.2326884687641787, 0.7596285680800611, 0.42707853717913047`],RGBColor[0.23984379062888905`, 0.7628837001518682, 0.4227860189043366],RGBColor[0.24715366801913152`, 0.7661133955136393, 0.4184127227273974],RGBColor[0.25461342821257693`, 0.7693184668521016, 0.4139573287646593],RGBColor[0.26221888931032444`, 0.7724987014653718, 0.40941845760197537`],RGBColor[0.26996266615436965`, 0.7756514499718966, 0.40479818880388174`],RGBColor[0.2778414184206625, 0.7787783143950607, 0.4000928616204019],RGBColor[0.2858480262057704, 0.7818756741313698, 0.3953036041833475],RGBColor[0.2939823104683633, 0.7849460858913447, 0.39042968981612],RGBColor[0.3022350954370879, 0.7879859567629643, 0.3854733248017807],RGBColor[0.310604418737109, 0.7909960529255057, 0.3804329133629437],RGBColor[0.3190871961690892, 0.7939771294497029, 0.37530572247730154`],RGBColor[0.32767842481037307`, 0.7969256090092027, 0.3700971086652205],RGBColor[0.3363770709617544, 0.7998440132840018, 0.3648008403071768],RGBColor[0.3451757700977104, 0.8027288102158318, 0.35942126510781447`],RGBColor[0.3540765354074894, 0.80558246359489, 0.353955095785225],RGBColor[0.36307288905435614`, 0.80840049484089, 0.348402577877059],RGBColor[0.3721639204075886, 0.8111846130302013, 0.3427637682343151],RGBColor[0.3813479050465821, 0.8139344946962004, 0.3370376142646351],RGBColor[0.3906182241219263, 0.8166487205595967, 0.33122632945843844`],RGBColor[0.3999768212623815, 0.8193276503392856, 0.32532725843019283`],RGBColor[0.40941597344689007`, 0.8219699071711365, 0.3193425718381236],RGBColor[0.418938095542911, 0.8245747957604673, 0.31327047427719934`],RGBColor[0.42853611297072447`, 0.8271430102244819, 0.30711089679069037`],RGBColor[0.43821148499483115`, 0.8296721588204928, 0.300864481689657],RGBColor[0.4479617886413208, 0.8321628664603256, 0.2945294594583805],RGBColor[0.4577813910404001, 0.8346148603683272, 0.2881066671370873],RGBColor[0.46767297073952857`, 0.8370283618478455, 0.2815945980811821],RGBColor[0.47762895833183777`, 0.8394011274151628, 0.27499396246121266`],RGBColor[0.4876519966897123, 0.8417343479104805, 0.2683017620352695],RGBColor[0.49773367432610355`, 0.8440258066445298, 0.2615237686948751],RGBColor[0.5078750219432514, 0.846277098641311, 0.2546561143755624],RGBColor[0.5180752494862431, 0.8484887693085125, 0.24769889915726157`],RGBColor[0.5283258321751821, 0.8506586693743721, 0.24065414755310766`],RGBColor[0.5386308275986138, 0.8527868884051059, 0.23351750058916276`],RGBColor[0.548980458002187, 0.8548743348948882, 0.22629568729131472`],RGBColor[0.5593788376902379, 0.8569210584030033, 0.21898246061294346`],RGBColor[0.5698174276705644, 0.8589260001543768, 0.21158392883280225`],RGBColor[0.5802945527724039, 0.8608896727270254, 0.2040979693981519],RGBColor[0.5908116383216248, 0.8628135638060236, 0.19652278866190287`],RGBColor[0.6013580840397512, 0.8646966714690255, 0.18886521134133236`],RGBColor[0.6119370906314997, 0.8665399511423325, 0.181121291501226],RGBColor[0.6225401569353771, 0.8683424453868425, 0.1732986740171394],RGBColor[0.6331682881172852, 0.8701070824773419, 0.16539484611862115`],RGBColor[0.6438129884912803, 0.8718329499665033, 0.1574190301339476],RGBColor[0.6544728151822515, 0.8735204919838195, 0.14937219444590882`],RGBColor[0.6651461868847677, 0.8751721361118794, 0.1412584133407773],RGBColor[0.6758257194413555, 0.8767870195906821, 0.13309191286806318`],RGBColor[0.6865115265930126, 0.8783679817337392, 0.1248770398795781],RGBColor[0.6971941434439218, 0.8799142035724945, 0.11663311070005217`],RGBColor[0.7078748175056854, 0.881428482695378, 0.10837392920767788`],RGBColor[0.7185439321356931, 0.882912055742582, 0.10012807215776806`],RGBColor[0.7292039995694638, 0.8843656662392944, 0.09192536612477004],RGBColor[0.7398433313314975, 0.8857905919469689, 0.08380928749785493],RGBColor[0.7504625667873018, 0.8871902151408768, 0.07583905286235423],RGBColor[0.7610594940142685, 0.8885638719723012, 0.0680900722152191],RGBColor[0.7716234421245648, 0.8899148892894538, 0.06066995495444857],RGBColor[0.7821578637018753, 0.8912449579301015, 0.05371196587603636],RGBColor[0.7926533490577626, 0.8925554087193903, 0.047401117840516725`],RGBColor[0.8031101670648537, 0.8938489258169644, 0.04195435979844082],RGBColor[0.8135200000367577, 0.8951268587699225, 0.03764562984062277],RGBColor[0.8238814958821538, 0.8963915533178293, 0.034852646118004794`],RGBColor[0.8341953189844913, 0.8976453379379864, 0.033642700571396934`],RGBColor[0.8444522250765154, 0.8988895853765866, 0.034045498308027174`],RGBColor[0.8546555615949485, 0.9001259343493543, 0.036085858895558834`],RGBColor[0.8647960127489654, 0.9013577843655304, 0.03978910711211607],RGBColor[0.8748769440175387, 0.9025857610712735, 0.04495402504133204],RGBColor[0.8848859474442141, 0.903815286843018, 0.05123011464126393],RGBColor[0.8948269727230432, 0.905044655982066, 0.05838119124087505],RGBColor[0.9047025295790283, 0.9062782052563443, 0.06620564078749713],RGBColor[0.9145057699880221, 0.9075153058386868, 0.07453488746258119],RGBColor[0.9242397319919776, 0.9087596065041207, 0.0832510020840733],RGBColor[0.9338985888012371, 0.910011479184347, 0.09225150797640028],RGBColor[0.9434864385344409, 0.9112725641897205, 0.10147371721450577`],RGBColor[0.9529994532916154, 0.9125452328290099, 0.11085876909361342`]}; (* The viridis color palette *)
(*colorFunc=Blend[Reverse[viridisColors],(Abs[#]+1*^-16-mm[[1]])/(mm[[2]]-mm[[1]])]&; Trying out auto scaling, if uncommenting this, add the option ColorFunctionScaling\[Rule]False to the plotting function.*)
colorFunc=If[#<0,0,Blend[Reverse[viridisColors],(#-mm[[1]])/mm[[2]]]]&;

arr=Table[None,{i,(vMax-vMin)*sc+1},{j,(vMax-vMin)*sc+1}]; (* The array (or matrix) of signal values to plot *)
For[i=1, i<= Length[tempSig],i++,
arr=ReplacePart[arr,
{Floor[(tempSig[[i]][[2]]-vMin)*sc+1],
Floor[(tempSig[[i]][[1]]-vMin)*sc+1]}->If[tempSig[[i]][[3]]<0,Pink,tempSig[[i]][[3]]/scale]];
];


graphs=ArrayPlot[Reverse[arr],
(* Graph Exterior *)
PlotLabel->pl,
FrameLabel->{"West Deflector (V)","East Deflector (V)"},
FrameTicks->ft,

PlotLegends->BarLegend[{colorFunc,{0,mm[[2]]}},LegendLabel->Superscript["x 10",ToString[ me[[2]] ]]],
(*
PlotLegends\[Rule]Automatic,
*)
DataRange->ConstantArray[{vMin,vMax},2],
ColorFunction->colorFunc,
ColorFunctionScaling->False
];

graphs
]



End[];


EndPackage[];

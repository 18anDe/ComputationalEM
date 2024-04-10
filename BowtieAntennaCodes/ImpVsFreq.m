antennaLength = 0.5;
noOfSteps = 200;
betaLStart = 5;
betaLStop = 270;
FreqStart = betaLStart*c0/(360*antennaLength);
FreqStop = betaLStop*c0/(360*antennaLength);
step = (FreqStop-FreqStart)/(noOfSteps-1);
f = zeros;
for FF = 1:noOfSteps
    f(FF) = FreqStart+step*(FF-1);
    L = 360*antennaLength*f/c0;
end

Imp_90 = 0.5*[0.232478700002369 - 1164.76342276932i	0.373724536229685 - 914.364492779602i	0.549131128568967 - 749.926219045401i	0.759312667892187 - 633.265712157125i	1.00500635406616 - 545.894596559549i	1.28707562877116 - 477.771068974302i	1.60651395630988 - 422.970420795484i	1.96444916560760 - 377.770954256204i	2.36214836793733 - 339.716912033489i	2.80102346599585 - 307.122219375920i	3.28263727074299 - 278.790702333761i	3.80871024281937 - 253.850080803980i	4.38112787528154 - 231.649185261026i	5.00194873371999 - 211.691983115622i	5.67341316841673 - 193.593902457279i	6.39795271088132 - 177.052131322657i	7.17820016367554 - 161.824940761852i	8.01700038764759 - 147.716988701774i	8.91742178425794 - 134.568680692467i	9.88276846223574 - 122.248340093444i	10.9165930669466 - 110.646360440265i	12.0227102370843 - 99.6707801587658i	13.2052106360493 - 89.2438937778812i	14.4684754839670 - 79.2996292587587i	15.8171914899424 - 69.7814990918830i	17.2563660519303 - 60.6409864334125i	18.7913425524918 - 51.8362649624030i	20.4278155314798 - 43.3311776124240i	22.1718454600349 - 35.0944183084014i	24.0298727726450 - 27.0988746092010i	26.0087307337639 - 19.3210992602355i	28.1156566208156 - 11.7408861497621i	30.3583005943931 - 4.34093176594235i	32.7447314971786 + 2.89343252297224i	35.2834386736145 + 9.97444883685035i	37.9833287309492 + 16.9120908302436i	40.8537159675509 + 23.7141785885298i	43.9043049756077 + 30.3864561388694i	47.1451636827232 + 36.9326362962832i	50.5866848321457 + 43.3544168161913i	54.2395336181247 + 49.6514713282215i	58.1145788976851 + 55.8214182532907i	62.2228051031147 + 61.8597708402337i	66.5752016957000 + 67.7598715914980i	71.1826267517435 + 73.5128146778804i	76.0556410851580 + 79.1073604708326i	81.2043092242004 + 84.5298470503499i	86.6379636205310 + 89.7641044752587i	92.3649287350287 + 94.7913787238358i	98.3922021859434 + 99.5902735063786i	104.725091040337 + 104.136719580673i	111.366802665647 + 108.403982704507i	118.317991422413 + 112.362722843005i	125.576264952304 + 115.981118580954i	133.135656957697 + 119.225071693831i	140.986077202664 + 122.058507283030i	149.112753953120 + 124.443784511819i	157.495689093917 + 126.342231487780i	166.109151479988 + 127.714814911754i	174.921239330651 + 128.522950464543i	183.893547148010 + 128.729453317926i	192.980976077162 + 128.299619562252i	202.131728062082 + 127.202418877134i	211.287522774927 + 125.411766854404i	220.384071353505 + 122.907832767259i	229.351831911574 + 119.678326369045i	238.117058382812 + 115.719696907954i	246.603136822834 + 111.038170556580i	254.732182748167 + 105.650550468163i	262.426851017028 + 99.5847080082838i	269.612288308277 + 92.8797051424306i	276.218139944683 + 85.5855064438412i	282.180510165526 + 77.7622637033492i	287.443770116750 + 69.4791846693212i	291.962112121695 + 60.8130272143602i	295.700762426549 + 51.8462879718574i	298.636786529736 + 42.6651770181052i	300.759449190313 + 33.3574848761901i	302.070122213026 + 24.0104533771543i	302.581763714330 + 14.7087573871155i	302.318019554456 + 5.53269103386825i	301.312018412391 - 3.44336815060107i	299.604945016847 - 12.1521681276304i	297.244480881434 - 20.5344408724836i	294.283199119965 - 28.5394113127995i	290.776990929000 - 36.1249781750418i	286.783587980637 - 43.2576040639796i	282.361229262286 - 49.9119625204726i	277.567504672590 - 56.0703949484588i	272.458392426314 - 61.7222308329989i	267.087494083324 - 66.8630215603536i	261.505460391207 - 71.4937324903283i	255.759593311317 - 75.6199307744840i	249.893604468562 - 79.2509986726659i	243.947507503765 - 82.3993945067616i	237.957620979380 - 85.0799764048677i	231.956659125886 - 87.3093979334795i	225.973889370037 - 89.1055797288786i	220.035337866116 - 90.4872573351322i	214.164026839057 - 91.4736025639417i	208.380230202098 - 92.0839136883172i	202.701736460699 - 92.3373685218024i	197.144110249366 - 92.2528337670509i	191.720945908039 - 91.8487238012179i	186.444108266051 - 91.1429021757430i	181.323957266777 - 90.1526194405044i	176.369554254108 - 88.8944813724331i	171.588848682257 - 87.3844422318752i	166.988844736954 - 85.6378182382347i	162.575747903638 - 83.6693170160595i	158.355091920317 - 81.4930792916538i	154.331846839585 - 79.1227296055055i	150.510509122557 - 76.5714332410625i	146.895174819826 - 73.8519569546596i	143.489596979561 - 70.9767314267787i	140.297228475228 - 67.9579136463915i	137.321251476575 - 64.8074476935743i	134.564594805493 - 61.5371226079348i	132.029940428545 - 58.1586262282868i	129.719720343432 - 54.6835940688376i	127.636105118442 - 51.1236524647060i	125.780985341464 - 47.4904553793704i	124.155947226792 - 43.7957144221747i	122.762243611028 - 40.0512217774523i	121.600761541108 - 36.2688658987568i	120.671987614585 - 32.4606399715059i	119.975972172335 - 28.6386432928860i	119.512293364663 - 24.8150758559735i	119.280022012544 - 21.0022265513567i	119.277688066686 - 17.2124555092620i	119.503249330209 - 13.4581711931837i	119.954062959361 - 9.75180291713899i	120.626860095851 - 6.10576948863405i	121.517723820597 - 2.53244467477364i	122.622070459548 + 0.955879852198348i	123.934634126308 + 4.34703351135733i	125.449454262465 + 7.62900728538171i	127.159865843558 + 10.7899953451887i	129.058491864791 + 13.8184376700957i	131.137237712659 + 16.7030636885117i	133.387287071963 + 19.4329371505609i	135.799099115008 + 21.9975026226770i	138.362406871377 + 24.3866341516582i	141.066216879630 + 26.5906867677397i	143.898810470309 + 28.6005515688530i	146.847747312923 + 30.4077151378955i	149.899872164603 + 32.0043239796360i	153.041326068004 + 33.3832545143700i	156.257563540208 + 34.5381889253918i	159.533377550132 + 35.4636968251306i	162.852934273878 + 36.1553222840399i	166.199819719882 + 36.6096752670803i	169.557100303334 + 36.8245259621617i	172.907399299465 + 36.7988998880075i	176.232990800183 + 36.5331710679463i	179.515912327733 + 36.0291499905012i	182.738096621839 + 35.2901625919307i	185.881522324954 + 34.3211161376657i	188.928382368952 + 33.1285476960071i	191.861267856115 + 31.7206509308279i	194.663364181461 + 30.1072772234281i	197.318655128078 + 28.2999076858117i	199.812129755742 + 26.3115934487443i	202.129986171261 + 24.1568626767451i	204.259825787653 + 21.8515940345268i	206.190831507271 + 19.4128577394344i	207.913923440930 + 16.8587267986010i	209.421886315456 + 14.2080624521166i	210.709463612000 + 11.4802791257514i	211.773414673416 + 8.69509524528540i	212.612532451278 + 5.87227700136991i	213.227621140759 + 3.03138252395561i	213.621434571341 + 0.191513902648342i	213.798577776866 - 2.62891592093375i	213.765375560893 - 5.41239512507730i	213.529713019749 - 8.14249863296481i	213.100853826538 - 10.8040434413712i	212.489242581644 - 13.3832110831460i	211.706297694408 - 15.8676360572974i	210.764201097197 - 18.2464604020217i	209.675690649156 - 20.5103557923211i	208.453860419468 - 22.6515155599565i	207.111973214214 - 24.6636198248799i	205.663288793914 - 26.5417774794814i	204.120910283518 - 28.2824490837832i	202.497650357908 - 29.8833548296696i	200.805917937543 - 31.3433716446846i	199.057625382387 - 32.6624232667633i	197.264115547209 - 33.8413667688256i	195.436107566050 - 34.8818785843485i	193.583659867366 - 35.7863426162748i	191.716148676034 - 36.5577425317070i	189.842260121026 - 37.1995598776926i	187.969994022162 - 37.7156792171146i	186.106677458641 - 38.1103010905261i	184.258986308907 - 38.3878632671888i	182.432973079828 - 38.5529704594583i	180.634099498805 - 38.6103324389753i	178.867272513317 - 38.5647103080113i	177.136882518428 - 38.4208705404254i];
Imp_70 = 0.5*[0.226850930486187 - 1489.41504703649i	0.364537556738196 - 1170.46238962815i	0.535378395605046 - 961.243395186444i	0.739878482460596 - 813.017281740534i	0.978643280156323 - 702.183784717779i	1.25238095292081 - 615.924889278955i	1.56190502235855 - 546.678147615688i	1.90813741344066 - 489.693304114656i	2.29211189920431 - 441.836101378176i	2.71497795357504 - 400.954793741353i	3.17800502226664 - 365.522997424936i	3.68258722205881 - 334.427776622266i	4.23024847885924 - 306.838435739088i	4.82264811476949 - 282.122299549518i	5.46158689383053 - 259.788955974150i	6.14901353514513 - 239.452338353913i	6.88703170057272 - 220.804326136649i	7.67790746205839 - 203.595979417614i	8.52407725076289 - 187.623951305292i	9.42815628636313 - 172.720485601927i	10.3929474800034 - 158.745943648353i	11.4214507982089 - 145.583145543761i	12.5168730673706 - 133.133033028423i	13.6826381889078 - 121.311308695154i	14.9223977235861 - 110.045805787862i	16.2400417893448 - 99.2744112717533i	17.6397101999563 - 88.9434125938551i	19.1258037514009 - 79.0061723231243i	20.7029955384756 - 69.4220590633675i	22.3762421552371 - 60.1555805906241i	24.1507945987330 - 51.1756780472410i	26.0322086553693 - 42.4551495736793i	28.0263545023662 - 33.9701789065210i	30.1394252022359 - 25.6999498693648i	32.3779437051427 - 17.6263317940004i	34.7487679015051 - 9.73362406184932i	37.2590931843553 - 2.00835038856063i	39.9164518870151 + 5.56090464018493i	42.7287088559185 + 12.9836228090826i	45.7040523005429 + 20.2674481820279i	48.8509789324261 + 27.4182556970149i	52.1782722637239 + 34.4401958020406i	55.6949727840909 + 41.3357180001859i	59.4103385753004 + 48.1055758357050i	63.3337947597595 + 54.7488156501859i	67.4748700175983 + 61.2627513564723i	71.8431182551235 + 67.6429275121254i	76.4480233757503 + 73.8830731189945i	81.2988850069173 + 79.9750488300824i	86.4046829905587 + 85.9087906071973i	91.7739184723758 + 91.6722533413022i	97.4144295527926 + 97.2513585167862i	103.333179721127 + 102.629950661988i	109.536017719336 + 107.789768064896i	116.027408111042 + 112.710434019236i	122.810132705178 + 117.369475664804i	129.884964139421 + 121.742378243433i	137.250314398966 + 125.802683239219i	144.901862851262 + 129.522139318623i	152.832170518769 + 132.870915125645i	161.030289762753 + 135.817882693514i	169.481381247983 + 138.330979371095i	178.166352889599 + 140.377654590181i	187.061538283516 + 141.925405391868i	196.138434664560 + 142.942401290276i	205.363522438747 + 143.398194739018i	214.698189466335 + 143.264508220049i	224.098783174023 + 142.516082944769i	233.516811899488 + 141.131567621143i	242.899313327134 + 139.094419117429i	252.189402279271 + 136.393780695739i	261.327002472265 + 133.025298456909i	270.249757345616 + 128.991833441318i	278.894104194482 + 124.304026148807i	287.196484313730 + 118.980672625206i	295.094650647840 + 113.048877030271i	302.529024628241 + 106.543954741123i	309.444046576011 + 99.5090721593487i	315.789460226781 + 91.9946236850759i	321.521472282806 + 84.0573616633924i	326.603732687846 + 75.7593101384446i	331.008090339639 + 67.1665065516450i	334.715091500570 + 58.3476258066034i	337.714203124841 + 49.3725474169906i	340.003759307398 + 40.3109282204587i	341.590644609286 + 31.2308403560957i	342.489741769606 + 22.1975273330662i	342.723182196478 + 13.2723209380914i	342.319444950265 + 4.51174957730484i	341.312353439987 - 4.03314434189975i	339.740018911493 - 12.3172730405351i	337.643776514030 - 20.3017644086234i	335.067154017466 - 27.9540311941632i	332.054905947969 - 35.2476710132565i	328.652137841210 - 42.1622242857612i	324.903537211449 - 48.6828192715348i	320.852720284760 - 54.7997334101888i	316.541696953417 - 60.5078985735320i	312.010451011769 - 65.8063750821214i	307.296628616145 - 70.6978158568119i	302.435325036579 - 75.1879382609675i	297.458958016655 - 79.2850173532108i	292.397215259242 - 82.9994106479338i	287.277063518260 - 86.3431212271072i	282.122807308513 - 89.3294032503261i	276.956186170800 - 91.9724116039802i	271.796500596474 - 94.2868956073141i	266.660758000939 - 96.2879353176472i	261.563831444264 - 97.9907179970671i	256.518625060647 - 99.4103516583827i	251.536241330785 - 100.561712237990i	246.626146384798 - 101.459320788928i	241.796330445341 - 102.117247096771i	237.053461308292 - 102.549036248573i	232.403029416947 - 102.767654893343i	227.849483624770 - 102.785454191132i	223.396357173591 - 102.614146733175i	219.046383752507 - 102.264795009666i	214.801603761017 - 101.747809292036i	210.663461091186 - 101.072953073868i	206.632890879823 - 100.249354473324i	202.710398773276 - 99.2855222366703i	198.896132303770 - 98.1893651958920i	195.189945005152 - 96.9682142232568i	191.591453903908 - 95.6288458927347i	188.100091013978 - 94.1775072039389i	184.715149445461 - 92.6199408502638i	181.435824711420 - 90.9614106213083i	178.261251786345 - 89.2067266222214i	175.190538436652 - 87.3602700714965i	172.222795309418 - 85.4260175056664i	169.357163231668 - 83.4075642761656i	166.592838139769 - 81.3081472719399i	163.929094027473 - 79.1306668426882i	161.365304272230 - 76.8777079332432i	158.900961672810 - 74.5515604708887i	156.535697507094 - 72.1542390752885i	154.269299897047 - 69.6875021864796i	152.101731748264 - 67.1528707307504i	150.033148513853 - 64.5516464682755i	148.063916016468 - 61.8849301908954i	146.194628547614 - 59.1536399642090i	144.426127449518 - 56.3585296361295i	142.759520371172 - 53.5002078649212i	141.196201375941 - 50.5791579544434i	139.737872062531 - 47.5957588236424i	138.386563842850 - 44.5503074822028i	137.144661498278 - 41.4430434355086i	136.014928108041 - 38.2741755006151i	135.000531408044 - 35.0439115816119i	134.105071592767 - 31.7524920283477i	133.332610513678 - 28.4002272875645i	132.687702151134 - 24.9875406505033i	132.175424138116 - 21.5150170056870i	131.801409987752 - 17.9834586194438i	131.571881515425 - 14.3939490877870i	131.493680742221 - 10.7479267287352i	131.574300310018 - 7.04726880895891i	131.821911118588 - 3.29438811522051i	132.245385499519 + 0.507656521555504i	132.854313757168 + 4.35503407941248i	133.659011319641 + 8.24299790000283i	134.670513040010 + 12.1657233498180i	135.900550359215 + 16.1161251992437i	137.361506082329 + 20.0856537130426i	139.066340433318 + 24.0640690683186i	141.028480859555 + 28.0391946578374i	143.261666798261 + 31.9966511916231i	145.779739368111 + 35.9195753853757i	148.596364830957 + 39.7883295353859i	151.724679861727 + 43.5802115381835i	155.176846424792 + 47.2691790089944i	158.963504726491 + 50.8256061250144i	163.093114735187 + 54.2160976125581i	167.571180659086 + 57.4033907061789i	172.399359131665 + 60.3463825072043i	177.574461247830 + 63.0003262369878i	183.087371465357 + 65.3172443299104i	188.921922875123 + 67.2466076576021i	195.053788048149 + 68.7363265727478i	201.449466371477 + 69.7340888587258i	208.065470202174 + 70.1890601422054i	214.847829798347 + 70.0539325836403i	221.732046210830 + 69.2872677525822i	228.643616889060 + 67.8560316572089i	235.499235753668 + 65.7381687576814i	242.208724758116 + 62.9250151079601i	248.677687748540 + 59.4233183718535i	254.810794782001 + 55.2566246994203i	260.515516336223 + 50.4658176538943i	265.706046737677 + 45.1086559308154i	270.307101017035 + 39.2582506086355i	274.257253337092 + 33.0005372092654i	277.511515713338 + 26.4309146449917i	280.042931150689 + 19.6503214240895i	281.843063938947 + 12.7610811943464i	282.921392785501 + 5.86286413728174i	283.303728253861 - 0.950922738627301i	283.029865935200 - 7.59607558906226i	282.150739182292 - 13.9997183187194i];
Imp_50 = 0.5*[0.217612487620252 - 1897.33423768688i	0.349595376940164 - 1492.24723455407i	0.513258690130298 - 1226.76306374693i	0.709022307481076 - 1038.87520105186i	0.937389282105714 - 898.560334290124i	1.19894756290044 - 789.512089088746i	1.49437200859203 - 702.110103376176i	1.82442669973960 - 630.311618440332i	2.18996755644585 - 570.129705984469i	2.59194527038908 - 518.827359148923i	3.03140856060172 - 474.463132739313i	3.50950776318201 - 435.621549191753i	4.02749876581706 - 401.246176869780i	4.58674729859952 - 370.532485055202i	5.18873359311437 - 342.856908119830i	5.83505742212773 - 317.728604326317i	6.52744353239641 - 294.755867665595i	7.26774748308911 - 273.622250827068i	8.05796190202933 - 254.069274736072i	8.90022317136896 - 235.883698636257i	9.79681855331994 - 218.888007032842i	10.7501937651278 - 202.933204078130i	11.7629610104662 - 187.893288494368i	12.8379074717622 - 173.660969612077i	13.9780042644901 - 160.144311792330i	15.1864158500443 - 147.264081543290i	16.4665098982473 - 134.951632358361i	17.8218675836515 - 123.147205256448i	19.2562942913064 - 111.798553788083i	20.7738306973297 - 100.859824600487i	22.3787641770812 - 90.2906410344921i	24.0756404786715 - 80.0553493690930i	25.8692755814673 - 70.1223964198608i	27.7647676377580 - 60.4638140668728i	29.7675088702382 - 51.0547915233310i	31.8831972678687 - 41.8733201802930i	34.1178478873021 - 32.8998989791100i	36.4778035257019 - 24.1172906943774i	38.9697444825894 - 15.5103214189307i	41.6006970725342 - 7.06571705008214i	44.3780404860982 + 1.22802822743307i	47.3095115225976 + 9.38075551388007i	50.4032066340265 + 17.4007201155304i	53.6675806240912 + 25.2946485042011i	57.1114412390232 + 33.0677755298393i	60.7439387672316 + 40.7238637319532i	64.5745496327828 + 48.2652063195537i	68.6130528235797 + 55.6926151952479i	72.8694978400077 + 63.0053952809966i	77.3541626858389 + 70.2013063578978i	82.0775002535712 + 77.2765136594698i	87.0500712861268 + 84.2255285586944i	92.2824619328890 + 91.0411408671727i	97.7851837699094 + 97.7143445235197i	103.568554034355 + 104.234258793314i	109.642553747981 + 110.588047538066i	116.016661393842 + 116.760839639191i	122.699659889333 + 122.735654285430i	129.699414796425 + 128.493335544534i	137.022622060720 + 134.012501432662i	144.674524112887 + 139.269513549475i	152.658593940042 + 144.238474233829i	160.976187782163 + 148.891259071134i	169.626168468559 + 153.197593389237i	178.604503112398 + 157.125182036042i	187.903840944916 + 160.639902140778i	197.513079489657 + 163.706068602734i	207.416930015153 + 166.286781590880i	217.595496183802 + 168.344364228437i	228.023882906480 + 169.840896731478i	238.671855428104 + 170.738850438318i	249.503571357372 + 171.001821311907i	260.477410401639 + 170.595357586966i	271.545927614401 + 169.487870323940i	282.655955622384 + 167.651608896889i	293.748879198979 + 165.063676196891i	304.761101380946 + 161.707051041256i	315.624713899284 + 157.571578551306i	326.268376012551 + 152.654883818553i	336.618395133197 + 146.963160795309i	346.599990445011 + 140.511787762216i	356.138707824737 + 133.325723543691i	365.161941848496 + 125.439645212252i	373.600509683642 + 116.897798339750i	381.390213461807 + 107.753544492741i	388.473323384561 + 98.0686067721226i	394.799914109718 + 87.9120315247222i	400.328992239719 + 77.3589014104659i	405.029362787024 + 66.4888502151657i	408.880196586397 + 55.3844417024203i	411.871277553050 + 44.1294822674901i	414.002926923084 + 32.8073395369706i	415.285619512939 + 21.4993362576554i	415.739323046898 + 10.2832812990521i	415.392604461806 - 0.767811709146818i	414.281555963667 - 11.5867814156967i	412.448598126278 - 22.1131297423779i	409.941217605724 - 32.2935798708672i	406.810693612881 - 42.0823876502864i	403.110860934500 - 51.4414240283637i	398.896948961144 - 60.3400545583920i	394.224526821493 - 68.7548475874169i	389.148575201832 - 76.6691455073647i	383.722696450837 - 84.0725337848119i	377.998466641815 - 90.9602408457413i	372.024926696176 - 97.3324988041353i	365.848204591993 - 103.193891005699i	359.511257073270 - 108.552707874672i	353.053717015738 - 113.420327986072i	346.511831504792 - 117.810636928912i	339.918475518935 - 121.739492577364i	333.303226662455 - 125.224241965254i	326.692487442549 - 128.283292115365i	320.109642954238 - 130.935734908804i	313.575243368881 - 133.201024355257i	307.107202199950 - 135.098703383827i	300.721002855829 - 136.648176446844i	294.429907423876 - 137.868523742512i	288.245162926317 - 138.778352646443i	282.176201427888 - 139.395681932824i	276.230831352479 - 139.737854507724i	270.415418185787 - 139.821474622674i	264.735053413789 - 139.662365848458i	259.193711087162 - 139.275546436948i	253.794391825659 - 138.675219060345i	248.539254400434 - 137.874772275478i	243.429735272191 - 136.886791404609i	238.466656633467 - 135.723076845721i	233.650323617559 - 134.394668120149i	228.980611405915 - 132.911872231471i	224.457043000172 - 131.284295146494i	220.078858432621 - 129.520875417576i	215.845076176510 - 127.629919147087i	211.754547490891 - 125.619135651476i	207.806004397972 - 123.495673316592i	203.998101947759 - 121.266155249795i	200.329455377832 - 118.936714430548i	196.798672727387 - 116.513028141615i	193.404383415782 - 114.000351530000i	190.145263247804 - 111.403550202061i	187.020056261480 - 108.727131802654i	184.027593789946 - 105.975276565106i	181.166811066972 - 103.151866848713i	178.436761666316 - 100.260515704365i	175.836630028138 - 97.3045945277968i	173.365742291176 - 94.2872598747220i	171.023575617076 - 91.2114795234418i	168.809766163063 - 88.0800578788474i	166.724115830723 - 84.8956608177327i	164.766597891875 - 81.6608400792054i	162.937361567085 - 78.3780573060470i	161.236735608031 - 75.0497078433549i	159.665230911557 - 71.6781443997877i	158.223542170454 - 68.2657006743046i	156.912548543758 - 64.8147150473743i	155.733313307412 - 61.3275544303119i	154.687082424300 - 57.8066383593485i	153.775281950934 - 54.2544634122617i	152.999514176304 - 50.6736280145837i	152.361552366619 - 47.0668576892618i	151.863333967862 - 43.4370307879200i	151.506952096491 - 39.7872047231621i	151.294645127243 - 36.1206426992274i	151.228784166374 - 32.4408409124565i	151.311858178975 - 28.7515561628209i	151.546456520965 - 25.0568337829668i	151.935248610519 - 21.3610357512243i	152.480960461007 - 17.6688688095147i	153.186347788928 - 13.9854123557928i	154.054165407089 - 10.3161458232006i	155.087132616816 - 6.66697519463885i	156.287894324908 - 3.04425823196106i	157.658977633125 + 0.545172075984744i	159.202743682330 + 4.09398642199107i	160.921334581966 + 7.59434110462805i	162.816615320623 + 11.0378637117811i	164.890110637262 + 14.4156433921169i	167.142936937091 + 17.7182266241065i	169.575729463151 + 20.9356194661020i	172.188565085485 + 24.0572973331025i	174.980881245080 + 27.0722233908283i	177.951391789233 + 29.9688766796127i	181.098000657077 + 32.7352910724958i	184.417714615653 + 35.3591061266824i	187.906556503255 + 37.8276307977584i	191.559480700881 + 40.1279208444286i	195.370292815223 + 42.2468705517366i	199.331575806228 + 44.1713191372658i	203.434625014835 + 45.8881718748142i	207.669394726217 + 47.3845355728346i	212.024459023101 + 48.6478675839874i	216.486989724191 + 49.6661370048562i	221.042754146289 + 50.4279961639276i	225.676135258789 + 50.9229599090494i	230.370176502316 + 51.1415896161363i	235.106653110692 + 51.0756782770283i	239.866171204825 + 50.7184325182685i	244.628295224176 + 50.0646469890884i	249.371703441087 + 49.1108662710107i	254.074370389883 + 47.8555293365115i];
Imp_30 = 0.5*[0.202604274506366 - 2474.91942560533i	0.325412088594788 - 1947.95671047484i	0.477622591104917 - 1602.87605158638i	0.659578342785281 - 1358.89326997682i	0.871689606704135 - 1176.89318954127i	1.11443567491914 - 1035.63160482287i	1.38836642135633 - 922.575109721726i	1.69410408724482 - 829.851337331738i	2.03234530638418 - 752.266332406976i	2.40386337847080 - 686.254579863800i	2.80951079968082 - 629.287045482755i	3.25022206069855 - 579.519936626035i	3.72701672339293 - 535.577228381584i	4.24100278837758 - 496.411068959728i	4.79338036674209 - 461.209359619253i	5.38544567030902 - 429.332901795292i	6.01859533584849 - 400.271634489548i	6.69433109976304 - 373.613523377458i	7.41426484083306 - 349.022030784588i	8.18012400967141 - 326.219526906238i	8.99375746456685 - 304.974891617505i	9.85714173437389 - 285.094121994400i	10.7723877300152 - 266.413128731801i	11.7417479269665 - 248.792148898152i	12.7676240417639 - 232.111367518560i	13.8525752260604 - 216.267453872284i	14.9993268020157 - 201.170797488427i	16.2107795627582 - 186.743284777171i	17.4900196612460 - 172.916497331660i	18.8403291099662 - 159.630242018889i	20.2651969124504 - 146.831344310182i	21.7683308454099 - 134.472652113342i	23.3536699072440 - 122.512209204199i	25.0253974445662 - 110.912566296829i	26.7879549629961 - 99.6402046065523i	28.6460566215180 - 88.6650519965387i	30.6047044008844 - 77.9600758559153i	32.6692039254932 - 67.5009400232957i	34.8451809044373 - 57.2657155597480i	37.1385981405234 - 47.2346371470108i	39.5557730353713 - 37.3898984597736i	42.1033954935756 - 27.7154811240768i	44.7885460985175 - 18.1970128955196i	47.6187143958546 - 8.82165152300849i	50.6018170769412 + 0.422008554485787i	53.7462158022595 + 9.54400896828867i	57.0607343430145 + 18.5530814613132i	60.5546746458823 + 27.4566846457508i	64.2378313398361 + 36.2610245887001i	68.1205041031990 + 44.9710600674157i	72.2135071916658 + 53.5904931130905i	76.5281752919501 + 62.1217452812112i	81.0763647089233 + 70.5659199455884i	85.8704487146135 + 78.9227508108565i	90.9233056834316 + 87.1905367767331i	96.2482984081269 + 95.3660632694634i	101.859242734481 + 103.444510187264i	107.770363369886 + 111.419346694495i	113.996234413462 + 119.282213253318i	120.551701827023 + 127.022791513232i	127.451784723439 + 134.628663001923i	134.711552002123 + 142.085157990974i	142.345970525308 + 149.375196464408i	150.369720724478 + 156.479123814774i	158.796975282122 + 163.374544748910i	167.641136387071 + 170.036159919508i	176.914527059636 + 176.435611022079i	186.628032244655 + 182.541341514514i	196.790685847846 + 188.318481722881i	207.409200727211 + 193.728768869438i	218.487439940895 + 198.730514452428i	230.025829397360 + 203.278633346515i	242.020714555064 + 207.324750865060i	254.463667069727 + 210.817405671579i	267.340751356548 + 213.702367637839i	280.631765947816 + 215.923090254057i	294.309480241043 + 217.421316685474i	308.338893610661 + 218.137856680031i	322.676550633524 + 218.013547887361i	337.269952937264 + 216.990409392363i	352.057114332248 + 215.012987113825i	366.966310669353 + 212.029880030687i	381.916078366403 + 207.995423062677i	396.815514766420 + 202.871487246856i	411.564928447499 + 196.629341400211i	426.056877474464 + 189.251502948361i	440.177617886025 + 180.733490648625i	453.808963500922 + 171.085380475036i	466.830532161106 + 160.333060028268i	479.122324375627 + 148.519078408710i	490.567550352748 + 135.702998993347i	501.055593617844 + 121.961182608069i	510.484977142104 + 107.385957692395i	518.766184350655 + 92.0841704497343i	525.824185079034 + 76.1751486802783i	531.600526889547 + 59.7881540888272i	536.054875017635 + 43.0594349942699i	539.165917835039 + 26.1290203915044i	540.931595869778 + 9.13741393610992i	541.368656869029 - 7.77764926792999i	540.511582496883 - 24.4842294152708i	538.410969701188 - 40.8588598979827i	535.131478166309 - 56.7887744538908i	530.749472546007 - 72.1736176044647i	525.350493806420 - 86.9266098161718i	519.026688870059 - 100.975173971728i	511.874313794767 - 114.261059725538i	503.991405558486 - 126.740025247617i	495.475693967627 - 138.381150926776i	486.422800890434 - 149.165866972854i	476.924751097769 - 159.086777438876i	467.068799019152 - 168.146358352108i	456.936559598061 - 176.355598967625i	446.603419486906 - 183.732644180232i	436.138196930968 - 190.301484204476i	425.603014406947 - 196.090725837531i	415.053346779711 - 201.132468730559i	404.538208727735 - 205.461300595587i	394.100447789562 - 209.113417404903i	383.777113000433 - 212.125868443077i	373.599873223590 - 214.535921457682i	363.595463550513 - 216.380539947991i	353.786142271956 - 217.695962613222i	344.190144724904 - 218.517373923778i	334.822123693853 - 218.878654460162i	325.693568942400 - 218.812199888096i	316.813200869257 - 218.348798036853i	308.187335246970 - 217.517554384301i	299.820217554152 - 216.345857219428i	291.714326604336 - 214.859374770847i	283.870648060436 - 213.082077600090i	276.288919055609 - 211.036280522097i	268.967845567393 - 208.742699207426i	261.905294455221 - 206.220517427869i	255.098462208575 - 203.487461623592i	248.544022495431 - 200.559880096058i	242.238254573359 - 197.452824670607i	236.177154549541 - 194.180133132147i	230.356531366957 - 190.754511124420i	224.772089264754 - 187.187612525950i	219.419498320758 - 183.490117582039i	214.294454540299 - 179.671808289826i	209.392730812915 - 175.741640709513i	204.710219920781 - 171.707814015911i	200.242970651853 - 167.577836216048i	195.987217948433 - 163.358586545902i	191.939407908524 - 159.056374626348i	188.096218353345 - 154.676996509203i	184.454575579440 - 150.225787781544i	181.011667827535 - 145.707673923092i	177.764955922164 - 141.127218129401i	174.712181465263 - 136.488666824672i	171.851372902972 - 131.795993093629i	169.180849726711 - 127.052938263171i	166.699225016846 - 122.263051862305i	164.405406488913 - 117.429730183876i	162.298596158074 - 112.556253664325i	160.378288696368 - 107.645823288401i	158.644268519071 - 102.701596214690i	157.096605600475 - 97.7267208050270i	155.735649985281 - 92.7243712263761i	154.562024929200 - 87.6977817773514i	153.576618571068 - 82.6502810731310i	152.780574008492 - 77.5853262017976i	152.175277619859 - 72.5065369416419i	151.762345447204 - 67.4177301025516i	151.543607427351 - 62.3229540245504i	151.521089232903 - 57.2265232325830i	151.696991460626 - 52.1330532081902i	152.073665882948 - 47.0474951951502i	152.653588459522 - 41.9751709073053i	153.439328790831 - 36.9218069517654i	154.433515685900 - 31.8935687195458i	155.638798512483 - 26.8970934278653i	157.057804002353 - 21.9395219239588i	158.693088198073 - 17.0285287794215i	160.547083253141 - 12.1723501172051i	162.622038836602 - 7.37980852145058i	164.919957948622 - 2.66033428446601i	167.442527027528 + 1.97601785250822i	170.191040323753 + 6.51855740531725i	173.166318634455 + 10.9559533334552i	176.368622636167 + 15.2762368261835i	179.797561223294 + 19.4668128794678i	183.451995458357 + 23.5144819083732i	187.329938965366 + 27.4054726673318i	191.428455849080 + 31.1254877454368i	195.743557497017 + 34.6597628574783i	200.270099912592 + 37.9931410555565i	205.001683529608 + 41.1101628326325i	209.930557760135 + 43.9951728706539i	215.047532817392 + 46.6324438958096i	220.341901617149 + 49.0063177379524i	225.801374778005 + 51.1013632489707i	231.412031892946 + 52.9025502186694i	237.158292311451 + 54.3954378439267i	243.022908632122 + 55.5663756701373i	248.986985941161 + 56.4027142522493i	255.030029525379 + 56.8930221004384i	261.130023328353 + 57.0273048136872i];
Imp_10 = 0.5*[0.173346859291728 - 3539.09084971567i	0.278360161908856 - 2787.93128978774i	0.408453176193550 - 2296.49686039352i	0.563878935015114 - 1949.43053896907i	0.744940401606570 - 1690.87748111208i	0.951991410074681 - 1490.50209036714i	1.18543776775639 - 1330.40754724402i	1.44573852465623 - 1199.35286347567i	1.73340741601509 - 1089.92174424060i	2.04901448491649 - 997.023257233429i	2.39318789273453 - 917.046530340986i	2.76661592617584 - 847.359190698569i	3.17004921067140 - 785.996816923646i	3.60430314093837 - 731.463600603499i	4.07026054066720 - 682.600371524063i	4.56887456449680 - 638.494843842587i	5.10117185673484 - 598.419122401670i	5.66825598266327 - 561.785275093540i	6.27131114975316 - 528.113158187303i	6.91160623770671 - 497.006725274304i	7.59049915795340 - 468.136319911407i	8.30944156507009 - 441.225259947153i	9.06998394457060 - 416.039547098386i	9.87378110363981 - 392.379884127783i	10.7225980936780 - 370.075417656790i	11.6183165959827 - 348.978786566602i	12.5629418045396 - 328.962168882253i	13.5586098427420 - 309.914099928124i	14.6075957539035 - 291.736891793164i	15.7123221086984 - 274.344525665518i	16.8753682761620 - 257.660919048211i	18.0994804086089 - 241.618492435371i	19.3875821948113 - 226.156976917276i	20.7427864399820 - 211.222416938866i	22.1684075355793 - 196.766332152914i	23.6679748866375 - 182.745009773305i	25.2452473692405 - 169.118904609986i	26.9042288958659 - 155.852128472222i	28.6491851715827 - 142.912014164240i	30.4846617294589 - 130.268742094202i	32.4155033389171 - 117.895019744745i	34.4468748860982 - 105.765806037358i	36.5842838303948 - 93.8580740630996i	38.8336043460461 - 82.1506068213566i	41.2011032617921 - 70.6238215659159i	43.6934679148054 - 59.2596191468469i	46.3178360370536 - 48.0412553936029i	49.0818277924553 - 36.9532321361401i	51.9935800810608 - 25.9812059297519i	55.0617832213090 - 15.1119129517666i	58.2957201122174 - 4.33310889017126i	61.7053079630500 + 6.36647704430683i	65.3011426571249 + 16.9971745695296i	69.0945457872288 + 27.5683915746480i	73.0976143603998 + 38.0886202384948i	77.3232731170289 + 48.5654298750194i	81.7853293399777 + 59.0054461292041i	86.4985299398119 + 69.4143159206465i	91.4786204874528 + 79.7966573041718i	96.7424057197624 + 90.1559931852741i	102.307810859727 + 100.494667590771i	108.193942862695 + 110.813742952780i	114.421150413523 + 121.112876616649i	121.011081144747 + 131.390174534988i	127.986734109345 + 141.642019864752i	135.372505007453 + 151.862873952408i	143.194221016455 + 162.045046984237i	151.479161288216 + 172.178435415161i	160.256058234053 + 182.250223194028i	169.555073594384 + 192.244543810195i	179.407741963366 + 202.142100341813i	189.846872888826 + 211.919741049294i	200.906400879387 + 221.549988705079i	212.621170618913 + 231.000522879724i	225.026642424900 + 240.233615935752i	238.158500529443 + 249.205525660710i	252.052144184323 + 257.865850474191i	266.742039024889 + 266.156857169570i	282.260903774174 + 274.012796419504i	298.638705531313 + 281.359228011300i	315.901435993744 + 288.112386191870i	334.069641593430 + 294.178625747050i	353.156683442578 + 299.454001554803i	373.166709136438 + 303.824048183097i	394.092328981881 + 307.163841190299i	415.912005373421 + 309.338437231511i	438.587187100867 + 310.203804386048i	462.059251436905 + 309.608365062815i	486.246356507951 + 307.395278318767i	511.040354345139 + 303.405582478602i	536.303969296013 + 297.482297913746i	561.868503197935 + 289.475548828489i	587.532381279364 + 279.248697611619i	613.060891676706 + 266.685393307772i	638.187484663398 + 251.697318112197i	662.616971723497 + 234.232278846586i	686.030886910736 + 214.282147220018i	708.095135271188 + 191.890023599909i	728.469855596298 + 167.155910166652i	746.821179662279 + 140.240159055155i	762.834303785161 + 111.364033893206i	776.227040227681 + 80.8069021585910i	786.762832292774 + 48.8998557183386i	794.262142126572 + 16.0159101709282i	798.611184734971 - 17.4426882147333i	799.767191743813 - 51.0591772327995i	797.759721008078 - 84.4178993879078i	792.687933719936 - 117.119609914047i	784.714173453583 - 148.795311334467i	774.054535054638 - 179.117658486499i	760.967353418706 - 207.809288645584i	745.740646212795 - 234.647790528315i	728.679512029652 - 259.467373951479i	710.094341819798 - 282.157590080162i	690.290486677162 - 302.659650794837i	669.559782654295 - 320.960995981838i	648.174101117792 - 337.088768462113i	626.380897935294 - 351.102798251982i	604.400590468960 - 363.088595947371i	582.425500606469 - 373.150733243521i	560.620059436792 - 381.406866310677i	539.121964672897 - 387.982548181873i	518.044004074308 - 393.006886725730i	497.476296182122 - 396.609037551214i	477.488744835264 - 398.915475295854i	458.133549745133 - 400.047959174766i	439.447657693814 - 400.122095568343i	421.455075407180 - 399.246397908976i	404.168994964485 - 397.521748723712i	387.593705811786 - 395.041177571513i	371.726284676861 - 391.889879636870i	356.558066837752 - 388.145411389064i	342.075910248073 - 383.878010993671i	328.263268894735 - 379.151001480603i	315.101094271445 - 374.021243735380i	302.568584671658 - 368.539614093223i	290.643801676594 - 362.751487713062i	279.304172151948 - 356.697214098824i	268.526892577886 - 350.412575263137i	258.289250840604 - 343.929220249975i	248.568878861689 - 337.275072200053i	239.343947733398 - 330.474705995073i	230.593315424826 - 323.549695876987i	222.296635659803 - 316.518933410473i	214.434435257000 - 309.398916827839i	206.988166068329 - 302.204013236572i	199.940236647413 - 294.946695437738i	193.274027915040 - 287.637755243665i	186.973896349931 - 280.286495231004i	181.025167607135 - 272.900900847955i	175.414122938508 - 265.487794732571i	170.127980347590 - 258.052975008615i	165.154872041964 - 250.601339218097i	160.483819439410 - 243.136995433422i	156.104706730105 - 235.663361973946i	152.008253787330 - 228.183257034950i	148.185989046279 - 220.698979425390i	144.630222828585 - 213.212381505447i	141.334021473499 - 205.724935317350i	138.291182540878 - 198.237792813555i	135.496211272467 - 190.751841005374i	132.944298432921 - 183.267752782291i	130.631299598158 - 175.786034087406i	128.553715913355 - 168.307068076959i	126.708676304371 - 160.831156841535i	125.093921092924 - 153.358561222414i	123.707786935844 - 145.889539218402i	122.549192981144 - 138.424383445348i	121.617628107172 - 130.963458082160i	120.913139084925 - 123.507235712536i	120.436319476647 - 116.056334450460i	120.188299055538 - 108.611555718724i	120.170733500756 - 101.173923033022i	120.385794088415 - 93.7447221283419i	120.836157062254 - 86.3255427488229i	121.524992326440 - 78.9183224058906i	122.455951057133 - 71.5253923913014i	123.633151778386 - 64.1495263103938i	125.061164391466 - 56.7939913750575i	126.744991584293 - 49.4626026640242i	128.690046979505 - 42.1597805182937i	130.902129305630 - 34.8906111896259i	133.387391796502 - 27.6609107979305i	136.152305939982 - 20.4772925762429i	139.203618609746 - 13.3472372869974i	142.548301524846 - 6.27916657731782i	146.193491893810 + 0.717481099327364i	150.146423016435 + 7.63217253603209i	154.414343541823 + 14.4532005407054i	159.004424021150 + 21.1676043949740i	163.923649355724 + 27.7610912516890i	169.178695733514 + 34.2179600316533i	174.775790681534 + 40.5210297343556i	180.720554949620 + 46.6515744725403i	187.017825097633 + 52.5892679735055i	193.671455899354 + 58.3121407554504i	200.684102019719 + 63.7965536747612i	208.056978886194 + 69.0171920353002i	215.789603278142 + 73.9470849325540i	223.879514917133 + 78.5576549461862i];

a=figure;
plot(L, real(Imp_90),'b','LineWidth',1.4);hold on;
plot(L, real(Imp_70),'g','LineWidth',1.4);hold on;
plot(L, real(Imp_50),'r','LineWidth',1.4);hold on;
plot(L, real(Imp_30),'m','LineWidth',1.4);hold on;
plot(L, real(Imp_10),'k','LineWidth',1.4);hold on;
xlabel ('Antenna Length, Degrees','FontSize',14)
ylabel('Input  resistance, Ohm','FontSize',14)
xlim([0 300])
ylim([0 450])
legend('\theta = 90','\theta = 70','\theta = 50','\theta = 30','\theta = 10','Location','northeast','FontSize',14);

b=figure;
plot(L, imag(Imp_90),'b','LineWidth',1.4);hold on;
plot(L, imag(Imp_70),'g','LineWidth',1.4);hold on;
plot(L, imag(Imp_50),'r','LineWidth',1.4);hold on;
plot(L, imag(Imp_30),'m','LineWidth',1.4);hold on;
plot(L, imag(Imp_10),'k','LineWidth',1.4);hold on;
xlabel ('Antenna Length, Degrees','FontSize',14)
ylabel('Input  reactance, Ohm','FontSize',14)
xlim([0 300])
ylim([-250 200])
legend('\theta = 90','\theta = 70','\theta = 50','\theta = 30','\theta = 10','Location','northeast','FontSize',14);

ax = gca;
ax.FontSize = 14;           % Font size adjusted to 14
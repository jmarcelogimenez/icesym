#### ---- ####
# Archivo generado por SimulatorGUI
# CIMEC - Santa Fe - Argentina 
# Adecuado para levantar desde Interfaz Grafica 
# O para correr desde consola mediante $python main.py 
#### ---- ####


#--------- Inicializacion de Simulator

Simulator0 = dict()
Simulator0['dtheta_rpm'] = 1.0
Simulator0['filein_state'] = ''
Simulator0['calc_engine_data'] = 1
Simulator0['Courant'] = 0.8
Simulator0['heat_flow'] = 1.0
Simulator0['R_gas'] = 287.0
Simulator0['rpms'] = [2800, 3000, 3200]
Simulator0['filesave_spd'] = ''
Simulator0['filesave_state'] = 'QUB400'
Simulator0['ncycles'] = 10
Simulator0['folder_name'] = 'QUB400'
Simulator0['ga'] = 1.4
Simulator0['viscous_flow'] = 1.0
Simulator0['nsave'] = 10
Simulator0['ig_order'] = [0, 1]
Simulator0['get_state'] = 2
Simulator0['nappend'] = 10.0
Simulator0['engine_type'] = 0
Simulator0['nstroke'] = 2

Simulator = Simulator0


#--------- FIN Inicializacion de Simulator


#--------- Inicializacion de Cylinders

Cylinders = []

Cylinders0 = dict()
Cylinders0['crank_radius'] = 0.035
Cylinders0['type_ig'] = 0
Cylinders0['ndof'] = 3
Cylinders0['full_implicit'] = 0.0
Cylinders0['model_ht'] = 1
Cylinders0['factor_ht'] = 0.7
Cylinders0['piston_area'] = 0.005674502
Cylinders0['ownState'] = 1
Cylinders0['mass_C'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
Cylinders0['nnod'] = 3
Cylinders0['label'] = 'cyl'
Cylinders0['scavenge_type'] = 0
Cylinders0['twall'] = [473.0]
Cylinders0['state_ini'] = [[1.1769, 101330.0, 300.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0]]
Cylinders0['nve'] = 1
Cylinders0['head_chamber_area'] = 0.005674502
Cylinders0['type_temperature'] = 0
Cylinders0['rod_length'] = 0.125
Cylinders0['species_model'] = 0
Cylinders0['nvi'] = 1
Cylinders0['delta_ca'] = 0.0
Cylinders0['extras'] = 1
Cylinders0['histo'] = [0]
Cylinders0['Vol_clearance'] = 4.3411488649e-05
Cylinders0['Bore'] = 0.085
Cylinders0['scavenge'] = 1.0


#--------- Inicializacion de fuel

fuel0 = dict()
fuel0['y'] = 2.25
fuel0['hvap_fuel'] = 350000.0
fuel0['Q_fuel'] = 43000000.0

Cylinders0['fuel'] = fuel0


#--------- FIN Inicializacion de fuel


#--------- Inicializacion de combustion

combustion0 = dict()
combustion0['phi'] = 0.9
combustion0['a_wiebe'] = 5.25
combustion0['dtheta_comb'] = 1.0471975512
combustion0['combustion_model'] = 1
combustion0['m_wiebe'] = 1.1
combustion0['theta_ig_0'] = 6.1261056745

Cylinders0['combustion'] = combustion0


#--------- FIN Inicializacion de combustion


#--------- Inicializacion de injection

injection0 = dict()

Cylinders0['injection'] = injection0


#--------- FIN Inicializacion de injection

Cylinders0['position'] = (328,133)
Cylinders0['exhaust_valves'] = []
Cylinders0['intake_valves'] = []

Cylinders.append(Cylinders0)

Cylinders1 = dict()
Cylinders1['crank_radius'] = 0.035
Cylinders1['type_ig'] = 0
Cylinders1['ndof'] = 3
Cylinders1['full_implicit'] = 0.0
Cylinders1['model_ht'] = 1
Cylinders1['factor_ht'] = 0.6
Cylinders1['piston_area'] = 0.0056745
Cylinders1['ownState'] = 1
Cylinders1['mass_C'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
Cylinders1['nnod'] = 3
Cylinders1['label'] = 'ccase'
Cylinders1['twall'] = [323.0]
Cylinders1['state_ini'] = [[1.1769, 101330.0, 300.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0]]
Cylinders1['nve'] = 1
Cylinders1['head_chamber_area'] = 0.0056745
Cylinders1['type_temperature'] = 0
Cylinders1['rod_length'] = 0.125
Cylinders1['species_model'] = 0
Cylinders1['nvi'] = 1
Cylinders1['delta_ca'] = 0.0
Cylinders1['extras'] = 1
Cylinders1['histo'] = [0]
Cylinders1['Vol_clearance'] = 0.000882700269196
Cylinders1['Bore'] = 0.085
Cylinders1['scavenge'] = 0.0


#--------- Inicializacion de fuel

fuel0 = dict()
fuel0['y'] = 2.25
fuel0['hvap_fuel'] = 350000.0
fuel0['Q_fuel'] = 1e-10

Cylinders1['fuel'] = fuel0


#--------- FIN Inicializacion de fuel


#--------- Inicializacion de combustion

combustion0 = dict()
combustion0['phi'] = 0.0
combustion0['a_wiebe'] = 6.02
combustion0['dtheta_comb'] = 0.0174532925199
combustion0['combustion_model'] = 1
combustion0['m_wiebe'] = 1.64
combustion0['theta_ig_0'] = 3.12413936107

Cylinders1['combustion'] = combustion0


#--------- FIN Inicializacion de combustion


#--------- Inicializacion de injection

injection0 = dict()

Cylinders1['injection'] = injection0


#--------- FIN Inicializacion de injection

Cylinders1['position'] = (326,278)
Cylinders1['exhaust_valves'] = []
Cylinders1['intake_valves'] = []

Cylinders.append(Cylinders1)


#--------- FIN Inicializacion de Cylinders


#--------- Inicializacion de Valves

Valves = []

Valves0 = dict()
Valves0['angle_VC'] = 4.22369678983
Valves0['Lvmax'] = 0.010179
Valves0['Lv'] = [[0.0, 0.0], [118.0, 1.21079192672e-17], [119.0, 0.000249011043032], [120.0, 0.000529430658714], [121.0, 0.000824492166162], [122.0, 0.00112718803044], [123.0, 0.00143295238367], [124.0, 0.00173834956705], [125.0, 0.00204057783855], [126.0, 0.00233788510203], [127.0, 0.00263005643894], [128.0, 0.00291705913446], [129.0, 0.00319886241812], [130.0, 0.00347543739953], [131.0, 0.00374675700403], [132.0, 0.0040127959081], [133.0, 0.00427353047495], [134.0, 0.00452893869012], [135.0, 0.0047790000974], [136.0, 0.00502369573514], [137.0, 0.00526300807289], [138.0, 0.00549692094874], [139.0, 0.00572541950715], [140.0, 0.00594849013766], [141.0, 0.00616612041421], [142.0, 0.00637829903553], [143.0, 0.00658501576635], [144.0, 0.00678626137965], [145.0, 0.00698202759999], [146.0, 0.00717230704802], [147.0, 0.00735709318605], [148.0, 0.007536380265], [149.0, 0.0077101632725], [150.0, 0.00787843788239], [151.0, 0.00804120040552], [152.0, 0.008198447742], [153.0, 0.00835013606886], [154.0, 0.00849604508771], [155.0, 0.00863594840655], [156.0, 0.00876967028928], [157.0, 0.00889707942659], [158.0, 0.00901808379531], [159.0, 0.00913262634957], [160.0, 0.00924068135456], [161.0, 0.00934225122543], [162.0, 0.00943736376899], [163.0, 0.00952606975122], [164.0, 0.00960844073113], [165.0, 0.00968456711464], [166.0, 0.00975455639208], [167.0, 0.00981853152972], [168.0, 0.00987662949181], [169.0, 0.00992899987329], [170.0, 0.00997580362722], [171.0, 0.0100172118731], [172.0, 0.0100534047745], [173.0, 0.0100845704764], [174.0, 0.0101109040929], [175.0, 0.0101326067382], [176.0, 0.0101498845942], [177.0, 0.010162948008], [178.0, 0.0101720106141], [179.0, 0.0101772884765], [180.0, 0.0101789992445], [181.0, 0.0101772884765], [182.0, 0.0101720106141], [183.0, 0.010162948008], [184.0, 0.0101498845942], [185.0, 0.0101326067382], [186.0, 0.0101109040929], [187.0, 0.0100845704764], [188.0, 0.0100534047745], [189.0, 0.0100172118731], [190.0, 0.00997580362722], [191.0, 0.00992899987329], [192.0, 0.00987662949181], [193.0, 0.00981853152972], [194.0, 0.00975455639208], [195.0, 0.00968456711464], [196.0, 0.00960844073113], [197.0, 0.00952606975122], [198.0, 0.00943736376899], [199.0, 0.00934225122543], [200.0, 0.00924068135456], [201.0, 0.00913262634957], [202.0, 0.00901808379531], [203.0, 0.00889707942659], [204.0, 0.00876967028928], [205.0, 0.00863594840655], [206.0, 0.00849604508771], [207.0, 0.00835013606886], [208.0, 0.008198447742], [209.0, 0.00804120040552], [210.0, 0.00787843788239], [211.0, 0.0077101632725], [212.0, 0.007536380265], [213.0, 0.00735709318605], [214.0, 0.00717230704802], [215.0, 0.00698202759999], [216.0, 0.00678626137965], [217.0, 0.00658501576635], [218.0, 0.00637829903553], [219.0, 0.00616612041421], [220.0, 0.00594849013766], [221.0, 0.00572541950715], [222.0, 0.00549692094874], [223.0, 0.00526300807289], [224.0, 0.00502369573514], [225.0, 0.0047790000974], [226.0, 0.00452893869012], [227.0, 0.00427353047495], [228.0, 0.0040127959081], [229.0, 0.00374675700403], [230.0, 0.00347543739953], [231.0, 0.00319886241812], [232.0, 0.00291705913446], [233.0, 0.00263005643894], [234.0, 0.00233788510203], [235.0, 0.00204057783855], [236.0, 0.00173834956705], [237.0, 0.00143295238367], [238.0, 0.00112718803044], [239.0, 0.000824492166162], [240.0, 0.000529430658714], [241.0, 0.000249011043032], [242.0, 6.05395953302e-18], [360.0, 0.0]]
Valves0['label'] = 'tport_c'
Valves0['histo'] = 0
Valves0['Nval'] = 1
Valves0['Dv'] = 0.040716
Valves0['type_dat'] = 1
Valves0['Cd'] = [[0.0, 0.3], [0.0055, 0.5], [0.011, 0.9]]
Valves0['type'] = 0
Valves0['angle_V0'] = 2.05948851735
Valves0['valve_model'] = 1
Valves0['position'] = (281,131)
Valves0['tube'] = 1
Valves0['ncyl'] = 0
Valves0['typeVal'] = 'int'
Cylinders0['intake_valves'].append(Valves0)
Valves.append(Valves0)

Valves1 = dict()
Valves1['angle_VC'] = 4.60766922527
Valves1['Lvmax'] = 0.0102605
Valves1['Lv'] = [[0.0, 0.0], [96.0, 0.0], [97.0, 0.000162807260319], [98.0, 0.000342490962168], [99.0, 0.000531915278506], [100.0, 0.000728322388918], [101.0, 0.000930009428436], [102.0, 0.00113575514369], [103.0, 0.00134461081095], [104.0, 0.00155580193478], [105.0, 0.00176867478421], [106.0, 0.00198266436297], [107.0, 0.00219727383831], [108.0, 0.00241206061245], [109.0, 0.00262662649038], [110.0, 0.00284061049944], [111.0, 0.00305368349451], [112.0, 0.0032655440054], [113.0, 0.00347591497208], [114.0, 0.00368454112936], [115.0, 0.00389118687627], [116.0, 0.00409563451346], [117.0, 0.00429768276416], [118.0, 0.00449714551672], [119.0, 0.00469385183566], [120.0, 0.00488769291467], [121.0, 0.00507863497299], [122.0, 0.00526665094605], [123.0, 0.00545171514192], [124.0, 0.0056338032057], [125.0, 0.00581287749794], [126.0, 0.00598884032802], [127.0, 0.00616158587493], [128.0, 0.00633101822384], [129.0, 0.00649705082102], [130.0, 0.00665960597425], [131.0, 0.00681861439278], [132.0, 0.00697401476098], [133.0, 0.0071257533416], [134.0, 0.00727378360456], [135.0, 0.00741806587841], [136.0, 0.00755856702183], [137.0, 0.007695260113], [138.0, 0.00782812415505], [139.0, 0.00795714379614], [140.0, 0.00808230906272], [141.0, 0.00820361510506], [142.0, 0.00832106195394], [143.0, 0.00843465428774], [144.0, 0.00854440120932], [145.0, 0.00865031603189], [146.0, 0.00875241607354], [147.0, 0.00885072245969], [148.0, 0.00894525993337], [149.0, 0.00903605667265], [150.0, 0.00912314411499], [151.0, 0.00920655678831], [152.0, 0.00928633214824], [153.0, 0.0093625104215], [154.0, 0.00943513445513], [155.0, 0.00950424957116], [156.0, 0.00956990342675], [157.0, 0.00963214587938], [158.0, 0.00969102885703], [159.0, 0.00974660623303], [160.0, 0.00979893370548], [161.0, 0.00984806868101], [162.0, 0.0098940701627], [163.0, 0.00993699864196], [164.0, 0.00997691599427], [165.0, 0.0100138853785], [166.0, 0.0100479711395], [167.0, 0.0100792387145], [168.0, 0.010107754542], [169.0, 0.0101335859738], [170.0, 0.0101568011899], [171.0, 0.0101774691158], [172.0, 0.0101956593425], [173.0, 0.0102114420484], [174.0, 0.0102248879232], [175.0, 0.0102360680944], [176.0, 0.0102450540544], [177.0, 0.0102519175902], [178.0, 0.010256730714], [179.0, 0.0102595655947], [180.0, 0.0102604944913], [181.0, 0.0102595655947], [182.0, 0.010256730714], [183.0, 0.0102519175902], [184.0, 0.0102450540544], [185.0, 0.0102360680944], [186.0, 0.0102248879232], [187.0, 0.0102114420484], [188.0, 0.0101956593425], [189.0, 0.0101774691158], [190.0, 0.0101568011899], [191.0, 0.0101335859738], [192.0, 0.010107754542], [193.0, 0.0100792387145], [194.0, 0.0100479711395], [195.0, 0.0100138853785], [196.0, 0.00997691599427], [197.0, 0.00993699864196], [198.0, 0.0098940701627], [199.0, 0.00984806868101], [200.0, 0.00979893370548], [201.0, 0.00974660623303], [202.0, 0.00969102885703], [203.0, 0.00963214587938], [204.0, 0.00956990342675], [205.0, 0.00950424957116], [206.0, 0.00943513445513], [207.0, 0.0093625104215], [208.0, 0.00928633214824], [209.0, 0.00920655678831], [210.0, 0.00912314411499], [211.0, 0.00903605667265], [212.0, 0.00894525993337], [213.0, 0.00885072245969], [214.0, 0.00875241607354], [215.0, 0.00865031603189], [216.0, 0.00854440120932], [217.0, 0.00843465428774], [218.0, 0.00832106195394], [219.0, 0.00820361510506], [220.0, 0.00808230906272], [221.0, 0.00795714379614], [222.0, 0.00782812415505], [223.0, 0.007695260113], [224.0, 0.00755856702183], [225.0, 0.00741806587841], [226.0, 0.00727378360456], [227.0, 0.0071257533416], [228.0, 0.00697401476098], [229.0, 0.00681861439278], [230.0, 0.00665960597425], [231.0, 0.00649705082102], [232.0, 0.00633101822384], [233.0, 0.00616158587493], [234.0, 0.00598884032802], [235.0, 0.00581287749794], [236.0, 0.0056338032057], [237.0, 0.00545171514192], [238.0, 0.00526665094605], [239.0, 0.00507863497299], [240.0, 0.00488769291467], [241.0, 0.00469385183566], [242.0, 0.00449714551672], [243.0, 0.00429768276416], [244.0, 0.00409563451346], [245.0, 0.00389118687627], [246.0, 0.00368454112936], [247.0, 0.00347591497208], [248.0, 0.0032655440054], [249.0, 0.00305368349451], [250.0, 0.00284061049944], [251.0, 0.00262662649038], [252.0, 0.00241206061245], [253.0, 0.00219727383831], [254.0, 0.00198266436297], [255.0, 0.00176867478421], [256.0, 0.00155580193478], [257.0, 0.00134461081095], [258.0, 0.00113575514369], [259.0, 0.000930009428436], [260.0, 0.000728322388918], [261.0, 0.000531915278506], [262.0, 0.000342490962168], [263.0, 0.000162807260319], [264.0, 0.0], [360.0, 0.0]]
Valves1['label'] = 'eport'
Valves1['histo'] = 0
Valves1['Nval'] = 1
Valves1['Dv'] = 0.041041978
Valves1['type_dat'] = 1
Valves1['Cd'] = [[0.0, 0.6], [0.0102604944913, 0.85]]
Valves1['type'] = 1
Valves1['angle_V0'] = 1.67551608191
Valves1['valve_model'] = 1
Valves1['position'] = (372,128)
Valves1['tube'] = 2
Valves1['ncyl'] = 0
Valves1['typeVal'] = 'exh'
Cylinders0['exhaust_valves'].append(Valves1)
Valves.append(Valves1)

Valves2 = dict()
Valves2['angle_VC'] = 3.13286600733
Valves2['Lv'] = [[0.0, 0.00876068767886], [360.0, 0.00876068767886]]
Valves2['label'] = 'tport_cc'
Valves2['histo'] = 0
Valves2['Nval'] = 1
Valves2['Dv'] = 0.03504275
Valves2['type_dat'] = 1
Valves2['Cd'] = [[0.0, 0.7], [0.00876068767886, 0.7]]
Valves2['type'] = 1
Valves2['angle_V0'] = 3.15031929985
Valves2['valve_model'] = 1
Valves2['position'] = (363,282)
Valves2['tube'] = 1
Valves2['ncyl'] = 1
Valves2['typeVal'] = 'exh'
Cylinders1['exhaust_valves'].append(Valves2)
Valves.append(Valves2)

Valves3 = dict()
Valves3['angle_VC'] = 4.18879020479
Valves3['Lvmax'] = 0.008577155
Valves3['Lv'] = [[0.0, 0.0], [120.0, 0.0], [121.0, 0.000204523342945], [122.0, 0.000423023115357], [123.0, 0.000648899502709], [124.0, 0.000879462320691], [125.0, 0.00111297200067], [126.0, 0.00134812440448], [127.0, 0.00158386162212], [128.0, 0.00181928226412], [129.0, 0.00205359204193], [130.0, 0.00228607348957], [131.0, 0.00251609180239], [132.0, 0.0027433447691], [133.0, 0.00296772671377], [134.0, 0.00318913589088], [135.0, 0.00340747149774], [136.0, 0.00362263374598], [137.0, 0.00383452393208], [138.0, 0.00404304450705], [139.0, 0.00424809914514], [140.0, 0.0044495928115], [141.0, 0.00464743182879], [142.0, 0.00484152394266], [143.0, 0.00503177838603], [144.0, 0.00521810594212], [145.0, 0.00540041900624], [146.0, 0.00557863164627], [147.0, 0.00575265966174], [148.0, 0.00592242064154], [149.0, 0.00608783402022], [150.0, 0.00624882087061], [151.0, 0.00640525679035], [152.0, 0.00655692830706], [153.0, 0.0067036349948], [154.0, 0.00684520191898], [155.0, 0.00698147773753], [156.0, 0.00711233309876], [157.0, 0.00723765926355], [158.0, 0.00735736689891], [159.0, 0.00747138500327], [160.0, 0.00757965993342], [161.0, 0.00768215451036], [162.0, 0.00777884718615], [163.0, 0.00786973125802], [164.0, 0.00795481411881], [165.0, 0.00803411653514], [166.0, 0.00810767194652], [167.0, 0.00817552578012], [168.0, 0.00823773477659], [169.0, 0.00829436632396], [170.0, 0.00834549779653], [171.0, 0.00839121589701], [172.0, 0.00843161600003], [173.0, 0.008466801496], [174.0, 0.00849688313424], [175.0, 0.00852197836497], [176.0, 0.00854221067948], [177.0, 0.00855770894849], [178.0, 0.00856860675845], [179.0, 0.00857504174583], [180.0, 0.00857715492955], [181.0, 0.00857504174583], [182.0, 0.00856860675845], [183.0, 0.00855770894849], [184.0, 0.00854221067948], [185.0, 0.00852197836497], [186.0, 0.00849688313424], [187.0, 0.008466801496], [188.0, 0.00843161600003], [189.0, 0.00839121589701], [190.0, 0.00834549779653], [191.0, 0.00829436632396], [192.0, 0.00823773477659], [193.0, 0.00817552578012], [194.0, 0.00810767194652], [195.0, 0.00803411653514], [196.0, 0.00795481411881], [197.0, 0.00786973125802], [198.0, 0.00777884718615], [199.0, 0.00768215451036], [200.0, 0.00757965993342], [201.0, 0.00747138500327], [202.0, 0.00735736689891], [203.0, 0.00723765926355], [204.0, 0.00711233309876], [205.0, 0.00698147773753], [206.0, 0.00684520191898], [207.0, 0.0067036349948], [208.0, 0.00655692830706], [209.0, 0.00640525679035], [210.0, 0.00624882087061], [211.0, 0.00608783402022], [212.0, 0.00592242064154], [213.0, 0.00575265966174], [214.0, 0.00557863164627], [215.0, 0.00540041900624], [216.0, 0.00521810594212], [217.0, 0.00503177838603], [218.0, 0.00484152394266], [219.0, 0.00464743182879], [220.0, 0.0044495928115], [221.0, 0.00424809914514], [222.0, 0.00404304450705], [223.0, 0.00383452393208], [224.0, 0.00362263374598], [225.0, 0.00340747149774], [226.0, 0.00318913589088], [227.0, 0.00296772671377], [228.0, 0.0027433447691], [229.0, 0.00251609180239], [230.0, 0.00228607348957], [231.0, 0.00205359204193], [232.0, 0.00181928226412], [233.0, 0.00158386162212], [234.0, 0.00134812440448], [235.0, 0.00111297200067], [236.0, 0.000879462320691], [237.0, 0.000648899502709], [238.0, 0.000423023115357], [239.0, 0.000204523342945], [240.0, 0.0], [360.0, 0.0]]
Valves3['label'] = 'iport'
Valves3['histo'] = 0
Valves3['Nval'] = 1
Valves3['Dv'] = 0.03430862
Valves3['type_dat'] = 1
Valves3['Cd'] = [[0.0, 0.7], [0.00857715492955, 0.9]]
Valves3['type'] = 0
Valves3['angle_V0'] = 2.09439510239
Valves3['valve_model'] = 1
Valves3['position'] = (285,281)
Valves3['tube'] = 0
Valves3['ncyl'] = 1
Valves3['typeVal'] = 'int'
Cylinders1['intake_valves'].append(Valves3)
Valves.append(Valves3)


#--------- FIN Inicializacion de Valves


#--------- Inicializacion de Tubes

Tubes = []

Tubes0 = dict()
Tubes0['diameter'] = [0.034, 0.034007715, 0.03401543, 0.034023145, 0.03403086, 0.034038575, 0.03404629, 0.034054005, 0.03406172, 0.034069435, 0.03407715, 0.034084865, 0.03409258, 0.034100295, 0.03410801, 0.034115725, 0.03412344, 0.034131155, 0.03413887, 0.034146585, 0.0341543, 0.034162015, 0.03416973, 0.034177445, 0.03418516, 0.034192875, 0.03420059, 0.034208305, 0.03421602, 0.034223735, 0.03423145, 0.034239165, 0.03424688, 0.034254595, 0.03426231, 0.034270025, 0.03427774, 0.034285455, 0.03429317, 0.034300885, 0.0343086]
Tubes0['numNorm'] = 2
Tubes0['nnod'] = 41
Tubes0['twall'] = [323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0, 323.0]
Tubes0['ndof'] = 3
Tubes0['state_ini'] = [[1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0]]
Tubes0['label'] = 'tube0'
Tubes0['histo'] = [0, 40]
Tubes0['xnod'] = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.2]
Tubes0['posNorm'] = [0.0, 1.0]
Tubes0['longitud'] = 0.2
Tubes0['typeSave'] = 1
Tubes0['position'] = (199,282)
Tubes0['nleft'] = 0
Tubes0['tleft'] = 'atmosphere'
Tubes0['nright'] = 1
Tubes0['tright'] = 'cylinder'

Tubes.append(Tubes0)

Tubes1 = dict()
Tubes1['diameter'] = [0.0350428, 0.0353579777778, 0.0356731555556, 0.0359883333333, 0.0363035111111, 0.0366186888889, 0.0369338666667, 0.0372490444444, 0.0375642222222, 0.0378794, 0.0381945777778, 0.0385097555556, 0.0388249333333, 0.0391401111111, 0.0394552888889, 0.0397704666667, 0.0400856444444, 0.0404008222222, 0.040716]
Tubes1['numNorm'] = 2
Tubes1['nnod'] = 19
Tubes1['twall'] = [343.0, 343.0, 343.0, 343.0, 343.0, 343.0, 343.0, 343.0, 343.0, 343.0, 343.0, 343.0, 343.0, 343.0, 343.0, 343.0, 343.0, 343.0, 343.0]
Tubes1['ndof'] = 3
Tubes1['state_ini'] = [[1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0]]
Tubes1['label'] = 'tube1'
Tubes1['histo'] = [0, 18]
Tubes1['xnod'] = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09]
Tubes1['posNorm'] = [0.0, 1.0]
Tubes1['longitud'] = 0.09
Tubes1['typeSave'] = 1
Tubes1['position'] = (318,212)
Tubes1['nleft'] = 1
Tubes1['tleft'] = 'cylinder'
Tubes1['nright'] = 0
Tubes1['tright'] = 'cylinder'

Tubes.append(Tubes1)

Tubes2 = dict()
Tubes2['diameter'] = [0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038]
Tubes2['numNorm'] = 2
Tubes2['nnod'] = 41
Tubes2['twall'] = [673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0, 673.0]
Tubes2['ndof'] = 3
Tubes2['state_ini'] = [[1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0]]
Tubes2['label'] = 'tube2'
Tubes2['histo'] = [0, 40]
Tubes2['xnod'] = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.2]
Tubes2['posNorm'] = [0.0, 1.0]
Tubes2['longitud'] = 0.2
Tubes2['typeSave'] = 1
Tubes2['position'] = (453,123)
Tubes2['nleft'] = 0
Tubes2['tleft'] = 'cylinder'
Tubes2['nright'] = 0
Tubes2['tright'] = 'tank'

Tubes.append(Tubes2)

Tubes3 = dict()
Tubes3['diameter'] = [0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038, 0.038]
Tubes3['numNorm'] = 2
Tubes3['nnod'] = 41
Tubes3['twall'] = [500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0]
Tubes3['ndof'] = 3
Tubes3['state_ini'] = [[1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0]]
Tubes3['label'] = 'tube3'
Tubes3['histo'] = [0, 40]
Tubes3['xnod'] = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.2]
Tubes3['posNorm'] = [0.0, 1.0]
Tubes3['longitud'] = 0.2
Tubes3['typeSave'] = 1
Tubes3['position'] = (595,136)
Tubes3['nleft'] = 0
Tubes3['tleft'] = 'tank'
Tubes3['nright'] = 1
Tubes3['tright'] = 'atmosphere'

Tubes.append(Tubes3)


#--------- FIN Inicializacion de Tubes


#--------- Inicializacion de Tanks

Tanks = []

Tanks0 = dict()
Tanks0['Area_wall'] = 0.0
Tanks0['nnod'] = 3
Tanks0['ndof'] = 3
Tanks0['label'] = 'tank0'
Tanks0['Volume'] = 0.01
Tanks0['state_ini'] = [[1.1769, 101330.0, 300.0], [1.1769, 0.1, 101330.0], [1.1769, 0.1, 101330.0]]
Tanks0['extras'] = 1
Tanks0['mass'] = 0.01769
Tanks0['histo'] = [0]
Tanks0['Cd_ports'] = [1.0, 1.0]
Tanks0['T_wall'] = 323.0
Tanks0['h_film'] = 0.0
Tanks0['position'] = (530,125)
Tanks0['int2tube'] = [2]
Tanks0['exh2tube'] = [3]

Tanks.append(Tanks0)


#--------- FIN Inicializacion de Tanks


#--------- Inicializacion de Junctions

Junctions = []


#--------- FIN Inicializacion de Junctions


#--------- Inicializacion de Atmospheres

Atmospheres = []

Atmospheres0 = dict()
Atmospheres0['nnod'] = 1
Atmospheres0['ndof'] = 3
Atmospheres0['state_ini'] = [1.1842, 0.1, 101330.0]
Atmospheres0['position'] = (127,275)

Atmospheres.append(Atmospheres0)

Atmospheres1 = dict()
Atmospheres1['nnod'] = 1
Atmospheres1['ndof'] = 3
Atmospheres1['state_ini'] = [1.1842, 0.1, 101330.0]
Atmospheres1['position'] = (687,125)

Atmospheres.append(Atmospheres1)


#--------- FIN Inicializacion de Atmospheres



kargs = {'Simulator':Simulator, 'Cylinders':Cylinders, 'Junctions':Junctions, 'Tubes':Tubes, 'Tanks':Tanks, 'Atmospheres':Atmospheres}

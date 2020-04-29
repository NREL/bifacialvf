import numpy as np

FIXED_ENDTOEND_GTIFRONT = np.array(
        [8.453096667259906, 55.175918755432214, 71.09405040555485, 108.75520006707735,
         124.70295079568395, 139.24391651348796, 111.90914237188711,
         115.25779271536798, 58.595470724661084, 19.492853393450336])

FIXED_ENDTOEND_GTIBACK = np.array(
        [1.6817631740125067, 10.2644974763074, 13.603638026567387, 20.331474588185692,
         23.392792734252055, 26.07716346329336, 21.304988511276385, 
         21.407148380581642, 11.262412057494464, 3.8770395811777])


SKY_BETA160_C05_D1 = np.array(
    [[0.39344519, 0.39292569, 0.39240611, 0.39188644, 0.3913667 ,
      0.39084688, 0.39032699, 0.38980703, 0.389287  , 0.38876691,
      0.38824676, 0.38772654, 0.38720627, 0.38668595, 0.38616558,
      0.38564515, 0.38512468, 0.38460417, 0.38408362, 0.38356304,
      0.38304241, 0.38252176, 0.38200108, 0.38148037, 0.38095963,
      0.38043888, 0.37991811, 0.37939732, 0.37887652, 0.37835571,
      0.37783489, 0.37731407, 0.37679325, 0.37627242, 0.3757516 ,
      0.37523079, 0.37470998, 0.37418919, 0.37366841, 0.37314765,
      0.37262691, 0.37210619, 0.37158549, 0.37106483, 0.37054419,
      0.37002359, 0.36950302, 0.36898249, 0.368462  , 0.36794156,
      0.36742116, 0.36690081, 0.36638051, 0.36586027, 0.36534008,
      0.36481996, 0.36429989, 0.3637799 , 0.36325996, 0.3627401 ,
      0.36222032, 0.3617006 , 0.36118097, 0.36066142, 0.36014195,
      0.35962256, 0.35910327, 0.35858406, 0.35806495, 0.35754594,
      0.35702703, 0.35650821, 0.35598951, 0.3554709 , 0.35495241,
      0.35443403, 0.35391576, 0.35339762, 0.35287959, 0.35236168,
      0.3518439 , 0.35132624, 0.35080872, 0.35029133, 0.34977407,
      0.34925695, 0.34873997, 0.34822313, 0.34770644, 0.34718989,
      0.3466735 , 0.34615726, 0.34564117, 0.34512524, 0.34460947,
      0.34409386, 0.34357842, 0.34306315, 0.34254804, 0.34203311],
     [0.39344519, 0.39292569, 0.39240611, 0.39188644, 0.3913667 ,
      0.39084688, 0.39032699, 0.38980703, 0.389287  , 0.38876691,
      0.38824676, 0.38772654, 0.38720627, 0.38668595, 0.38616558,
      0.38564515, 0.38512468, 0.38460417, 0.38408362, 0.38356304,
      0.38304241, 0.38252176, 0.38200108, 0.38148037, 0.38095963,
      0.38043888, 0.37991811, 0.37939732, 0.37887652, 0.37835571,
      0.37783489, 0.37731407, 0.37679325, 0.37627242, 0.3757516 ,
      0.37523079, 0.37470998, 0.37418919, 0.37366841, 0.37314765,
      0.37262691, 0.37210619, 0.37158549, 0.37106483, 0.37054419,
      0.37002359, 0.36950302, 0.36898249, 0.368462  , 0.36794156,
      0.36742116, 0.36690081, 0.36638051, 0.36586027, 0.36534008,
      0.36481996, 0.36429989, 0.3637799 , 0.36325996, 0.3627401 ,
      0.36222032, 0.3617006 , 0.36118097, 0.36066142, 0.36014195,
      0.35962256, 0.35910327, 0.35858406, 0.35806495, 0.35754594,
      0.35702703, 0.35650821, 0.35598951, 0.3554709 , 0.35495241,
      0.35443403, 0.35391576, 0.35339762, 0.35287959, 0.35236168,
      0.3518439 , 0.35132624, 0.35080872, 0.35029133, 0.34977407,
      0.34925695, 0.34873997, 0.34822313, 0.34770644, 0.34718989,
      0.3466735 , 0.34615726, 0.34564117, 0.34512524, 0.34460947,
      0.34409386, 0.34357842, 0.34306315, 0.34254804, 0.34203311]])

SKY_BETA20_C05_D1 = np.array(
    [[0.49191249, 0.47907243, 0.46639084, 0.4539596 , 0.44186595,
      0.43019095, 0.41900828, 0.40838328, 0.39837241, 0.38902304,
      0.38037348, 0.37245339, 0.36528422, 0.35888002, 0.35324812,
      0.34839005, 0.34430235, 0.34097738, 0.3384041 , 0.33656874,
      0.33545546, 0.3350468 , 0.3353242 , 0.33626834, 0.33785938,
      0.34007726, 0.34290179, 0.34631281, 0.35029021, 0.35460871,
      0.35899677, 0.36385087, 0.36914875, 0.374868  , 0.38098599,
      0.3874798 , 0.39432619, 0.40150151, 0.4089817 , 0.41674221,
      0.42475803, 0.43300368, 0.44145321, 0.45008021, 0.45885791,
      0.46775917, 0.4767566 , 0.4858226 , 0.4949295 , 0.50404957,
      0.51315521, 0.52221897, 0.53121369, 0.54011256, 0.5488892 ,
      0.55751777, 0.56597297, 0.57423011, 0.58226515, 0.59005469,
      0.59757597, 0.60480687, 0.61172586, 0.61831194, 0.62454461,
      0.63040379, 0.63586971, 0.64092289, 0.645544  , 0.64971378,
      0.65341298, 0.65662231, 0.65932232, 0.66149343, 0.66311588,
      0.66420495, 0.66492791, 0.66502809, 0.66448647, 0.66328438,
      0.6614038 , 0.65882772, 0.65554061, 0.65152894, 0.64678181,
      0.64129165, 0.63505501, 0.62807333, 0.62035387, 0.61191049,
      0.60276447, 0.59294524, 0.58249086, 0.57144843, 0.55987415,
      0.54783307, 0.53539856, 0.52265139, 0.5096785 , 0.49657147],
     [0.49191249, 0.47907243, 0.46639084, 0.4539596 , 0.44186595,
      0.43019095, 0.41900828, 0.40838328, 0.39837241, 0.38902304,
      0.38037348, 0.37245339, 0.36528422, 0.35888002, 0.35324812,
      0.34839005, 0.34430235, 0.34097738, 0.3384041 , 0.33656874,
      0.33545546, 0.3350468 , 0.3353242 , 0.33626834, 0.33785938,
      0.34007726, 0.34290179, 0.34631281, 0.35029021, 0.35460871,
      0.35899677, 0.36385087, 0.36914875, 0.374868  , 0.38098599,
      0.3874798 , 0.39432619, 0.40150151, 0.4089817 , 0.41674221,
      0.42475803, 0.43300368, 0.44145321, 0.45008021, 0.45885791,
      0.46775917, 0.4767566 , 0.4858226 , 0.4949295 , 0.50404957,
      0.51315521, 0.52221897, 0.53121369, 0.54011256, 0.5488892 ,
      0.55751777, 0.56597297, 0.57423011, 0.58226515, 0.59005469,
      0.59757597, 0.60480687, 0.61172586, 0.61831194, 0.62454461,
      0.63040379, 0.63586971, 0.64092289, 0.645544  , 0.64971378,
      0.65341298, 0.65662231, 0.65932232, 0.66149343, 0.66311588,
      0.66420495, 0.66492791, 0.66502809, 0.66448647, 0.66328438,
      0.6614038 , 0.65882772, 0.65554061, 0.65152894, 0.64678181,
      0.64129165, 0.63505501, 0.62807333, 0.62035387, 0.61191049,
      0.60276447, 0.59294524, 0.58249086, 0.57144843, 0.55987415,
      0.54783307, 0.53539856, 0.52265139, 0.5096785 , 0.49657147]])

SKY_BETA20_C0_D1 = np.array(
    [[0.02721499, 0.0283665 , 0.02958612, 0.03087902, 0.03225084,
      0.03370777, 0.03525655, 0.03690459, 0.03866   , 0.04053169,
      0.04252945, 0.04466404, 0.04694732, 0.04939235, 0.05201353,
      0.05482677, 0.05784964, 0.06110152, 0.06460389, 0.06838047,
      0.07245748, 0.07686393, 0.08163185, 0.08679662, 0.09239721,
      0.0984765 , 0.10508156, 0.11226383, 0.12007933, 0.12858868,
      0.13785703, 0.14795376, 0.15895183, 0.17092689, 0.18395573,
      0.19811434, 0.21347513, 0.23010344, 0.24805335, 0.26736258,
      0.28804702, 0.31009493, 0.3334614 , 0.35806381, 0.38377894,
      0.41044252, 0.4378518 , 0.46577152, 0.49394289, 0.52209507,
      0.54995789, 0.57727461, 0.60381305, 0.62937442, 0.65379876,
      0.67696726, 0.6988015 , 0.71926031, 0.73833509, 0.75604423,
      0.77242734, 0.78753973, 0.80144751, 0.81422338, 0.82594327,
      0.83668366, 0.84651965, 0.85552351, 0.86376381, 0.87130483,
      0.87820626, 0.88452317, 0.89030601, 0.89560085, 0.90044957,
      0.90489014, 0.90895687, 0.91268077, 0.91608977, 0.91920901,
      0.92206108, 0.92466627, 0.92704278, 0.9292069 , 0.93117321,
      0.93295473, 0.93456304, 0.93600842, 0.9373    , 0.93844578,
      0.93945279, 0.94032712, 0.941074  , 0.94169784, 0.94220232,
      0.94259038, 0.94286425, 0.94302554, 0.94307515, 0.94301338],
     [0.02721499, 0.0283665 , 0.02958612, 0.03087902, 0.03225084,
      0.03370777, 0.03525655, 0.03690459, 0.03866   , 0.04053169,
      0.04252945, 0.04466404, 0.04694732, 0.04939235, 0.05201353,
      0.05482677, 0.05784964, 0.06110152, 0.06460389, 0.06838047,
      0.07245748, 0.07686393, 0.08163185, 0.08679662, 0.09239721,
      0.0984765 , 0.10508156, 0.11226383, 0.12007933, 0.12858868,
      0.13785703, 0.14795376, 0.15895183, 0.17092689, 0.18395573,
      0.19811434, 0.21347513, 0.23010344, 0.24805335, 0.26736258,
      0.28804702, 0.31009493, 0.3334614 , 0.35806381, 0.38377894,
      0.41044252, 0.4378518 , 0.46577152, 0.49394289, 0.52209507,
      0.54995789, 0.57727461, 0.60381305, 0.62937442, 0.65379876,
      0.67696726, 0.6988015 , 0.71926031, 0.73833509, 0.75604423,
      0.77242734, 0.78753973, 0.80144751, 0.81422338, 0.82594327,
      0.83668366, 0.84651965, 0.85552351, 0.86376381, 0.87130483,
      0.87820626, 0.88452317, 0.89030601, 0.89560085, 0.90044957,
      0.90489014, 0.90895687, 0.91268077, 0.91608977, 0.91920901,
      0.92206108, 0.92466627, 0.92704278, 0.9292069 , 0.93117321,
      0.93295473, 0.93456304, 0.93600842, 0.9373    , 0.93844578,
      0.93945279, 0.94032712, 0.941074  , 0.94169784, 0.94220232,
      0.94259038, 0.94286425, 0.94302554, 0.94307515, 0.94301338]])

SKY_BETA160_C0_D1 = np.array(
    [[0.00384724, 0.00384054, 0.00383385, 0.00382718, 0.00382052,
      0.00381388, 0.00380726, 0.00380064, 0.00379405, 0.00378746,
      0.0037809 , 0.00377434, 0.0037678 , 0.00376128, 0.00375477,
      0.00374828, 0.0037418 , 0.00373533, 0.00372888, 0.00372244,
      0.00371602, 0.00370961, 0.00370321, 0.00369683, 0.00369047,
      0.00368412, 0.00367778, 0.00367145, 0.00366514, 0.00365885,
      0.00365256, 0.0036463 , 0.00364004, 0.0036338 , 0.00362757,
      0.00362136, 0.00361516, 0.00360898, 0.0036028 , 0.00359665,
      0.0035905 , 0.00358437, 0.00357825, 0.00357215, 0.00356606,
      0.00355998, 0.00355391, 0.00354786, 0.00354183, 0.0035358 ,
      0.00352979, 0.00352379, 0.00351781, 0.00351183, 0.00350588,
      0.00349993, 0.003494  , 0.00348808, 0.00348217, 0.00347628,
      0.0034704 , 0.00346453, 0.00345867, 0.00345283, 0.003447  ,
      0.00344118, 0.00343538, 0.00342958, 0.0034238 , 0.00341804,
      0.00341228, 0.00340654, 0.00340081, 0.0033951 , 0.00338939,
      0.0033837 , 0.00337802, 0.00337235, 0.0033667 , 0.00336105,
      0.00335542, 0.0033498 , 0.0033442 , 0.0033386 , 0.00333302,
      0.00332745, 0.00332189, 0.00331634, 0.00331081, 0.00330529,
      0.00329978, 0.00329428, 0.00328879, 0.00328332, 0.00327785,
      0.0032724 , 0.00326696, 0.00326153, 0.00325612, 0.00325071],
     [0.00384724, 0.00384054, 0.00383385, 0.00382718, 0.00382052,
      0.00381388, 0.00380726, 0.00380064, 0.00379405, 0.00378746,
      0.0037809 , 0.00377434, 0.0037678 , 0.00376128, 0.00375477,
      0.00374828, 0.0037418 , 0.00373533, 0.00372888, 0.00372244,
      0.00371602, 0.00370961, 0.00370321, 0.00369683, 0.00369047,
      0.00368412, 0.00367778, 0.00367145, 0.00366514, 0.00365885,
      0.00365256, 0.0036463 , 0.00364004, 0.0036338 , 0.00362757,
      0.00362136, 0.00361516, 0.00360898, 0.0036028 , 0.00359665,
      0.0035905 , 0.00358437, 0.00357825, 0.00357215, 0.00356606,
      0.00355998, 0.00355391, 0.00354786, 0.00354183, 0.0035358 ,
      0.00352979, 0.00352379, 0.00351781, 0.00351183, 0.00350588,
      0.00349993, 0.003494  , 0.00348808, 0.00348217, 0.00347628,
      0.0034704 , 0.00346453, 0.00345867, 0.00345283, 0.003447  ,
      0.00344118, 0.00343538, 0.00342958, 0.0034238 , 0.00341804,
      0.00341228, 0.00340654, 0.00340081, 0.0033951 , 0.00338939,
      0.0033837 , 0.00337802, 0.00337235, 0.0033667 , 0.00336105,
      0.00335542, 0.0033498 , 0.0033442 , 0.0033386 , 0.00333302,
      0.00332745, 0.00332189, 0.00331634, 0.00331081, 0.00330529,
      0.00329978, 0.00329428, 0.00328879, 0.00328332, 0.00327785,
      0.0032724 , 0.00326696, 0.00326153, 0.00325612, 0.00325071]])

SKY_BETA160_C1_D1 = np.array(
    [[0.3114026 , 0.31120886, 0.31101505, 0.31082117, 0.31062722,
      0.3104332 , 0.3102391 , 0.31004494, 0.30985071, 0.30965641,
      0.30946205, 0.30926761, 0.30907311, 0.30887854, 0.3086839 ,
      0.3084892 , 0.30829443, 0.30809959, 0.30790469, 0.30770973,
      0.3075147 , 0.3073196 , 0.30712444, 0.30692922, 0.30673394,
      0.30653859, 0.30634318, 0.30614771, 0.30595217, 0.30575658,
      0.30556092, 0.30536521, 0.30516943, 0.30497359, 0.3047777 ,
      0.30458174, 0.30438573, 0.30418966, 0.30399353, 0.30379734,
      0.30360109, 0.30340479, 0.30320844, 0.30301202, 0.30281555,
      0.30261903, 0.30242244, 0.30222581, 0.30202912, 0.30183238,
      0.30163558, 0.30143873, 0.30124183, 0.30104487, 0.30084786,
      0.3006508 , 0.30045369, 0.30025653, 0.30005932, 0.29986206,
      0.29966474, 0.29946738, 0.29926997, 0.29907251, 0.298875  ,
      0.29867745, 0.29847984, 0.29828219, 0.29808449, 0.29788675,
      0.29768895, 0.29749112, 0.29729323, 0.29709531, 0.29689733,
      0.29669932, 0.29650126, 0.29630315, 0.296105  , 0.29590681,
      0.29570858, 0.2955103 , 0.29531198, 0.29511362, 0.29491522,
      0.29471678, 0.2945183 , 0.29431978, 0.29412122, 0.29392262,
      0.29372398, 0.2935253 , 0.29332658, 0.29312783, 0.29292903,
      0.29273021, 0.29253134, 0.29233244, 0.2921335 , 0.29193452],
     [0.3114026 , 0.31120886, 0.31101505, 0.31082117, 0.31062722,
      0.3104332 , 0.3102391 , 0.31004494, 0.30985071, 0.30965641,
      0.30946205, 0.30926761, 0.30907311, 0.30887854, 0.3086839 ,
      0.3084892 , 0.30829443, 0.30809959, 0.30790469, 0.30770973,
      0.3075147 , 0.3073196 , 0.30712444, 0.30692922, 0.30673394,
      0.30653859, 0.30634318, 0.30614771, 0.30595217, 0.30575658,
      0.30556092, 0.30536521, 0.30516943, 0.30497359, 0.3047777 ,
      0.30458174, 0.30438573, 0.30418966, 0.30399353, 0.30379734,
      0.30360109, 0.30340479, 0.30320844, 0.30301202, 0.30281555,
      0.30261903, 0.30242244, 0.30222581, 0.30202912, 0.30183238,
      0.30163558, 0.30143873, 0.30124183, 0.30104487, 0.30084786,
      0.3006508 , 0.30045369, 0.30025653, 0.30005932, 0.29986206,
      0.29966474, 0.29946738, 0.29926997, 0.29907251, 0.298875  ,
      0.29867745, 0.29847984, 0.29828219, 0.29808449, 0.29788675,
      0.29768895, 0.29749112, 0.29729323, 0.29709531, 0.29689733,
      0.29669932, 0.29650126, 0.29630315, 0.296105  , 0.29590681,
      0.29570858, 0.2955103 , 0.29531198, 0.29511362, 0.29491522,
      0.29471678, 0.2945183 , 0.29431978, 0.29412122, 0.29392262,
      0.29372398, 0.2935253 , 0.29332658, 0.29312783, 0.29292903,
      0.29273021, 0.29253134, 0.29233244, 0.2921335 , 0.29193452]])

SKY_BETA20_C1_D1 = np.array(
    [[0.48359335, 0.48106107, 0.47853666, 0.47603262, 0.47356126,
      0.4711346 , 0.46876431, 0.46646169, 0.46423755, 0.46210222,
      0.46006544, 0.45813638, 0.45632354, 0.45463477, 0.45307722,
      0.45165731, 0.45038074, 0.44925246, 0.44827667, 0.44745684,
      0.44679565, 0.44629509, 0.4459564 , 0.44578011, 0.44576606,
      0.44591339, 0.44622062, 0.44668562, 0.44730565, 0.4480774 ,
      0.44899698, 0.45006001, 0.45126157, 0.4525963 , 0.45405837,
      0.45564154, 0.45733918, 0.45914431, 0.4610496 , 0.46304743,
      0.46512988, 0.4672888 , 0.46951581, 0.47180231, 0.47413957,
      0.47651868, 0.47893064, 0.48136632, 0.48381656, 0.48627214,
      0.48872384, 0.49116241, 0.49357868, 0.49596352, 0.49830788,
      0.50060282, 0.50283954, 0.50500939, 0.50708952, 0.50899944,
      0.51081437, 0.51252644, 0.51412804, 0.51561185, 0.51697089,
      0.51819852, 0.51928846, 0.52023486, 0.52103225, 0.52167564,
      0.5221605 , 0.52248279, 0.522639  , 0.52262617, 0.52244188,
      0.52208434, 0.52155235, 0.52084535, 0.51996342, 0.51890734,
      0.51767854, 0.51627917, 0.51471209, 0.51298085, 0.51108975,
      0.50904377, 0.50684862, 0.50451069, 0.50203706, 0.49943547,
      0.49671429, 0.49388246, 0.49094953, 0.48792553, 0.48482096,
      0.48164675, 0.47841417, 0.47513479, 0.47182042, 0.468483  ],
     [0.48359335, 0.48106107, 0.47853666, 0.47603262, 0.47356126,
      0.4711346 , 0.46876431, 0.46646169, 0.46423755, 0.46210222,
      0.46006544, 0.45813638, 0.45632354, 0.45463477, 0.45307722,
      0.45165731, 0.45038074, 0.44925246, 0.44827667, 0.44745684,
      0.44679565, 0.44629509, 0.4459564 , 0.44578011, 0.44576606,
      0.44591339, 0.44622062, 0.44668562, 0.44730565, 0.4480774 ,
      0.44899698, 0.45006001, 0.45126157, 0.4525963 , 0.45405837,
      0.45564154, 0.45733918, 0.45914431, 0.4610496 , 0.46304743,
      0.46512988, 0.4672888 , 0.46951581, 0.47180231, 0.47413957,
      0.47651868, 0.47893064, 0.48136632, 0.48381656, 0.48627214,
      0.48872384, 0.49116241, 0.49357868, 0.49596352, 0.49830788,
      0.50060282, 0.50283954, 0.50500939, 0.50708952, 0.50899944,
      0.51081437, 0.51252644, 0.51412804, 0.51561185, 0.51697089,
      0.51819852, 0.51928846, 0.52023486, 0.52103225, 0.52167564,
      0.5221605 , 0.52248279, 0.522639  , 0.52262617, 0.52244188,
      0.52208434, 0.52155235, 0.52084535, 0.51996342, 0.51890734,
      0.51767854, 0.51627917, 0.51471209, 0.51298085, 0.51108975,
      0.50904377, 0.50684862, 0.50451069, 0.50203706, 0.49943547,
      0.49671429, 0.49388246, 0.49094953, 0.48792553, 0.48482096,
      0.48164675, 0.47841417, 0.47513479, 0.47182042, 0.468483  ]])

SKY_BETA20_C1_D0 = np.array(
    [[0.09025658, 0.09055715, 0.09085192, 0.09114066, 0.09142313,
      0.09169909, 0.09196829, 0.09223048, 0.09248542, 0.09273283,
      0.09297247, 0.09320406, 0.09342733, 0.09364202, 0.09384785,
      0.09404453, 0.0942318 , 0.09440935, 0.09457692, 0.09473421,
      0.09488092, 0.09501677, 0.09514146, 0.0952547 , 0.09535618,
      0.09544562, 0.09552272, 0.09558717, 0.09563868, 0.09567695,
      0.09570169, 0.09571261, 0.09570941, 0.0956918 , 0.09565949,
      0.09561221, 0.09554967, 0.09547159, 0.09537771, 0.09526775,
      0.09514146, 0.09499859, 0.09483889, 0.09466211, 0.09446802,
      0.09425641, 0.09402706, 0.09377976, 0.09351431, 0.09323055,
      0.09292828, 0.09260736, 0.09226763, 0.09190897, 0.09153124,
      0.09113435, 0.0907182 , 0.09028272, 0.08982784, 0.08935352,
      0.08885973, 0.08834646, 0.08781372, 0.08726153, 0.08668994,
      0.086099  , 0.0854888 , 0.08485944, 0.08421102, 0.0835437 ,
      0.08285763, 0.08215298, 0.08142994, 0.08068875, 0.07992961,
      0.0791528 , 0.07835859, 0.07754725, 0.07671911, 0.07587448,
      0.07501371, 0.07413717, 0.07324522, 0.07233826, 0.07141669,
      0.07048094, 0.06953143, 0.06856862, 0.06759297, 0.06660493,
      0.065605  , 0.06459365, 0.06357138, 0.06253868, 0.06149607,
      0.06044405, 0.05938314, 0.05831385, 0.05723669, 0.05615218],
     [0.09025658, 0.09055715, 0.09085192, 0.09114066, 0.09142313,
      0.09169909, 0.09196829, 0.09223048, 0.09248542, 0.09273283,
      0.09297247, 0.09320406, 0.09342733, 0.09364202, 0.09384785,
      0.09404453, 0.0942318 , 0.09440935, 0.09457692, 0.09473421,
      0.09488092, 0.09501677, 0.09514146, 0.0952547 , 0.09535618,
      0.09544562, 0.09552272, 0.09558717, 0.09563868, 0.09567695,
      0.09570169, 0.09571261, 0.09570941, 0.0956918 , 0.09565949,
      0.09561221, 0.09554967, 0.09547159, 0.09537771, 0.09526775,
      0.09514146, 0.09499859, 0.09483889, 0.09466211, 0.09446802,
      0.09425641, 0.09402706, 0.09377976, 0.09351431, 0.09323055,
      0.09292828, 0.09260736, 0.09226763, 0.09190897, 0.09153124,
      0.09113435, 0.0907182 , 0.09028272, 0.08982784, 0.08935352,
      0.08885973, 0.08834646, 0.08781372, 0.08726153, 0.08668994,
      0.086099  , 0.0854888 , 0.08485944, 0.08421102, 0.0835437 ,
      0.08285763, 0.08215298, 0.08142994, 0.08068875, 0.07992961,
      0.0791528 , 0.07835859, 0.07754725, 0.07671911, 0.07587448,
      0.07501371, 0.07413717, 0.07324522, 0.07233826, 0.07141669,
      0.07048094, 0.06953143, 0.06856862, 0.06759297, 0.06660493,
      0.065605  , 0.06459365, 0.06357138, 0.06253868, 0.06149607,
      0.06044405, 0.05938314, 0.05831385, 0.05723669, 0.05615218]])

SKY_BETA160_C1_D0 = np.array(
    [[0.0005987 , 0.00179573, 0.00299166, 0.00418578, 0.00537735,
      0.00656566, 0.00774999, 0.00892963, 0.0101039 , 0.01127208,
      0.01243352, 0.01358753, 0.01473347, 0.0158707 , 0.01699857,
      0.01811649, 0.01922386, 0.02032009, 0.02140464, 0.02247695,
      0.0235365 , 0.02458279, 0.02561534, 0.02663368, 0.02763738,
      0.02862601, 0.02959919, 0.03055652, 0.03149767, 0.0324223 ,
      0.0333301 , 0.03422078, 0.0350941 , 0.0359498 , 0.03678767,
      0.0376075 , 0.03840914, 0.03919242, 0.03995722, 0.04070342,
      0.04143094, 0.0421397 , 0.04282965, 0.04350077, 0.04415304,
      0.04478647, 0.04540109, 0.04599692, 0.04657403, 0.0471325 ,
      0.0476724 , 0.04819385, 0.04869696, 0.04918185, 0.04964868,
      0.05009759, 0.05052876, 0.05094236, 0.05133857, 0.0517176 ,
      0.05207964, 0.05242492, 0.05275366, 0.05306608, 0.05336243,
      0.05364293, 0.05390785, 0.05415743, 0.05439194, 0.05461162,
      0.05481675, 0.0550076 , 0.05518443, 0.05534752, 0.05549714,
      0.05563356, 0.05575708, 0.05586795, 0.05596646, 0.05605289,
      0.05612752, 0.05619061, 0.05624245, 0.05628331, 0.05631347,
      0.05633319, 0.05634274, 0.05634239, 0.05633241, 0.05631306,
      0.05628459, 0.05624726, 0.05620133, 0.05614705, 0.05608466,
      0.05601441, 0.05593653, 0.05585127, 0.05575886, 0.05565952],
     [0.0005987 , 0.00179573, 0.00299166, 0.00418578, 0.00537735,
      0.00656566, 0.00774999, 0.00892963, 0.0101039 , 0.01127208,
      0.01243352, 0.01358753, 0.01473347, 0.0158707 , 0.01699857,
      0.01811649, 0.01922386, 0.02032009, 0.02140464, 0.02247695,
      0.0235365 , 0.02458279, 0.02561534, 0.02663368, 0.02763738,
      0.02862601, 0.02959919, 0.03055652, 0.03149767, 0.0324223 ,
      0.0333301 , 0.03422078, 0.0350941 , 0.0359498 , 0.03678767,
      0.0376075 , 0.03840914, 0.03919242, 0.03995722, 0.04070342,
      0.04143094, 0.0421397 , 0.04282965, 0.04350077, 0.04415304,
      0.04478647, 0.04540109, 0.04599692, 0.04657403, 0.0471325 ,
      0.0476724 , 0.04819385, 0.04869696, 0.04918185, 0.04964868,
      0.05009759, 0.05052876, 0.05094236, 0.05133857, 0.0517176 ,
      0.05207964, 0.05242492, 0.05275366, 0.05306608, 0.05336243,
      0.05364293, 0.05390785, 0.05415743, 0.05439194, 0.05461162,
      0.05481675, 0.0550076 , 0.05518443, 0.05534752, 0.05549714,
      0.05563356, 0.05575708, 0.05586795, 0.05596646, 0.05605289,
      0.05612752, 0.05619061, 0.05624245, 0.05628331, 0.05631347,
      0.05633319, 0.05634274, 0.05634239, 0.05633241, 0.05631306,
      0.05628459, 0.05624726, 0.05620133, 0.05614705, 0.05608466,
      0.05601441, 0.05593653, 0.05585127, 0.05575886, 0.05565952]])

SKY_BETA160_C05_D0 = np.array(
    [[0.0019083 , 0.00572098, 0.00952188, 0.01330324, 0.0170574 ,
      0.02077689, 0.02445443, 0.028083  , 0.03165587, 0.03516663,
      0.03860922, 0.04197798, 0.04526763, 0.04847333, 0.05159066,
      0.05461567, 0.05754483, 0.06037508, 0.0631038 , 0.06572882,
      0.06824839, 0.0706612 , 0.07296634, 0.07516327, 0.07725184,
      0.07923225, 0.08110501, 0.08287096, 0.08453121, 0.08608712,
      0.08754031, 0.08889258, 0.09014596, 0.09130261, 0.09236487,
      0.09333518, 0.09421611, 0.0950103 , 0.09572047, 0.0963494 ,
      0.0968999 , 0.0973748 , 0.09777697, 0.09810924, 0.09837448,
      0.09857549, 0.09871507, 0.09879598, 0.09882094, 0.09879261,
      0.09871361, 0.09858648, 0.09841371, 0.09819775, 0.09794093,
      0.09764555, 0.09731382, 0.09694789, 0.09654982, 0.09612162,
      0.09566521, 0.09518243, 0.09467507, 0.09414484, 0.09359337,
      0.09302222, 0.09243291, 0.09182686, 0.09120544, 0.09056996,
      0.08992166, 0.08926174, 0.08859131, 0.08791146, 0.08722319,
      0.08652749, 0.08582525, 0.08511736, 0.08440462, 0.08368783,
      0.08296772, 0.08224497, 0.08152024, 0.08079415, 0.08006728,
      0.07934017, 0.07861333, 0.07788726, 0.07716238, 0.07643914,
      0.07571792, 0.0749991 , 0.07428301, 0.07356998, 0.0728603 ,
      0.07215426, 0.07145211, 0.07075408, 0.0700604 , 0.06937127],
     [0.0019083 , 0.00572098, 0.00952188, 0.01330324, 0.0170574 ,
      0.02077689, 0.02445443, 0.028083  , 0.03165587, 0.03516663,
      0.03860922, 0.04197798, 0.04526763, 0.04847333, 0.05159066,
      0.05461567, 0.05754483, 0.06037508, 0.0631038 , 0.06572882,
      0.06824839, 0.0706612 , 0.07296634, 0.07516327, 0.07725184,
      0.07923225, 0.08110501, 0.08287096, 0.08453121, 0.08608712,
      0.08754031, 0.08889258, 0.09014596, 0.09130261, 0.09236487,
      0.09333518, 0.09421611, 0.0950103 , 0.09572047, 0.0963494 ,
      0.0968999 , 0.0973748 , 0.09777697, 0.09810924, 0.09837448,
      0.09857549, 0.09871507, 0.09879598, 0.09882094, 0.09879261,
      0.09871361, 0.09858648, 0.09841371, 0.09819775, 0.09794093,
      0.09764555, 0.09731382, 0.09694789, 0.09654982, 0.09612162,
      0.09566521, 0.09518243, 0.09467507, 0.09414484, 0.09359337,
      0.09302222, 0.09243291, 0.09182686, 0.09120544, 0.09056996,
      0.08992166, 0.08926174, 0.08859131, 0.08791146, 0.08722319,
      0.08652749, 0.08582525, 0.08511736, 0.08440462, 0.08368783,
      0.08296772, 0.08224497, 0.08152024, 0.08079415, 0.08006728,
      0.07934017, 0.07861333, 0.07788726, 0.07716238, 0.07643914,
      0.07571792, 0.0749991 , 0.07428301, 0.07356998, 0.0728603 ,
      0.07215426, 0.07145211, 0.07075408, 0.0700604 , 0.06937127]])

SKY_BETA20_C05_D0 = np.array(
    [[0.09228642, 0.09322975, 0.09418129, 0.09514089, 0.09610841,
      0.09708366, 0.09806647, 0.09905662, 0.10005388, 0.10105801,
      0.10206874, 0.10308576, 0.10410876, 0.1051374 , 0.1061713 ,
      0.10721007, 0.10825326, 0.10930043, 0.11035108, 0.11140467,
      0.11246065, 0.11351841, 0.11457731, 0.11563666, 0.11669575,
      0.1177538 , 0.11880999, 0.11986347, 0.12091331, 0.12195855,
      0.12299818, 0.12403112, 0.12505623, 0.12607232, 0.12707815,
      0.1280724 , 0.1290537 , 0.13002059, 0.13097158, 0.13190509,
      0.13281946, 0.133713  , 0.13458391, 0.13543033, 0.13625035,
      0.13704196, 0.13780311, 0.13853166, 0.1392254 , 0.13988208,
      0.14049937, 0.14107488, 0.14160618, 0.14209076, 0.14241963,
      0.14264508, 0.14281147, 0.14291603, 0.14295597, 0.14292848,
      0.14283078, 0.14266007, 0.1424136 , 0.14208863, 0.14168249,
      0.1411926 , 0.14061643, 0.13995158, 0.13919576, 0.13834685,
      0.13740287, 0.13636205, 0.1352228 , 0.1339838 , 0.13264396,
      0.1312025 , 0.1296589 , 0.12801301, 0.12626498, 0.12441538,
      0.12246511, 0.12041551, 0.11826831, 0.1160257 , 0.11369029,
      0.11126511, 0.10875367, 0.10615991, 0.10348821, 0.10074338,
      0.09793063, 0.09505559, 0.09212424, 0.08914292, 0.08611827,
      0.08305723, 0.07996696, 0.07685483, 0.07372836, 0.07059518],
     [0.09228642, 0.09322975, 0.09418129, 0.09514089, 0.09610841,
      0.09708366, 0.09806647, 0.09905662, 0.10005388, 0.10105801,
      0.10206874, 0.10308576, 0.10410876, 0.1051374 , 0.1061713 ,
      0.10721007, 0.10825326, 0.10930043, 0.11035108, 0.11140467,
      0.11246065, 0.11351841, 0.11457731, 0.11563666, 0.11669575,
      0.1177538 , 0.11880999, 0.11986347, 0.12091331, 0.12195855,
      0.12299818, 0.12403112, 0.12505623, 0.12607232, 0.12707815,
      0.1280724 , 0.1290537 , 0.13002059, 0.13097158, 0.13190509,
      0.13281946, 0.133713  , 0.13458391, 0.13543033, 0.13625035,
      0.13704196, 0.13780311, 0.13853166, 0.1392254 , 0.13988208,
      0.14049937, 0.14107488, 0.14160618, 0.14209076, 0.14241963,
      0.14264508, 0.14281147, 0.14291603, 0.14295597, 0.14292848,
      0.14283078, 0.14266007, 0.1424136 , 0.14208863, 0.14168249,
      0.1411926 , 0.14061643, 0.13995158, 0.13919576, 0.13834685,
      0.13740287, 0.13636205, 0.1352228 , 0.1339838 , 0.13264396,
      0.1312025 , 0.1296589 , 0.12801301, 0.12626498, 0.12441538,
      0.12246511, 0.12041551, 0.11826831, 0.1160257 , 0.11369029,
      0.11126511, 0.10875367, 0.10615991, 0.10348821, 0.10074338,
      0.09793063, 0.09505559, 0.09212424, 0.08914292, 0.08611827,
      0.08305723, 0.07996696, 0.07685483, 0.07372836, 0.07059518]])
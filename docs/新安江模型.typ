#import "@local/modern-cug-report:0.1.2": *
// #import "@preview/modern-cug-report:0.1.1": *
#counter(heading).update(2)
#let delta(x) = $Delta #x$

#show: (doc) => template(doc, 
  header: "水文模型原理")


#figure(
  image("images/三水源新安江模型结构及计算层次.jpg", width: 120%),
  caption: [新安江模型结构及计算层次]
) <fig_1>


超渗产流：雨强超过土壤下渗能力而产生的径流(气候干燥区)\

蓄满产流：土壤含水量饱和后而产生的径流(气候湿润区)\

== 1 蒸散发计算
三层蒸散发模型：
按土壤垂向分布的不均匀性将土层分为三层。

$
W M = U M+L M+D M \
W   = W U + W L +W D \
E   = E U + E L + E D \ 
E P = K C · E M 
$
W : 总张力水蓄量，mm （WU、WL、WD分别为上层、下层、深层）\
E : 总的蒸散发量，mm \
EP ：蒸散发能力\
三层蒸发模式按照先上层，后下层的次序，分以下几种情况：\
① 若 P+ W U ≥ EP ，则 EU =EP,EL = 0 ,ED=0
#figure(
  image("images/蒸发.jpg", width: 80%),
  caption: [Es/Ep与θ的关系，Es：土壤蒸发量，Ep：流域蒸发能力]
) <fig_2>

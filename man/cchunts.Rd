\name{cchunts}
\alias{cchunts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Human forager hunting returns data}
\description{
  Available data sets in package cchunts.
}
\usage{
data(Alvard)
data(Beckerman)
data(Bird_Bird_Codding)
data(Coad)
data(Duda)
data(Ellen)
data(Fernandez_Llamazares)
data(Franzen)
data(Gallois)
data(Gueze)
data(Headland)
data(Healey_Nen_PNG)
data(Healey)
data(Hill_Kintigh)
data(Koster)
data(Kramer_Greaves)
data(Lupo_Schmitt)
data(Marks)
data(Napitupulu)
data(Nielsen)
data(Pacheco)
data(Pangau_Adam)
data(Ready)
data(Reyes_Garcia)
data(Sillitoe)
data(Siren)
data(Trumble_Gurven)
data(Van_Vliet_et_al_Baego)
data(Van_Vliet_et_al_Djoutou)
data(Van_Vliet_et_al_Gabon)
data(Van_Vliet_et_al_Ingolo)
data(Van_Vliet_et_al_Ngombe)
data(Van_Vliet_et_al_Ovan)
data(Van_Vliet_et_al_Phalanga)
data(Van_Vliet_et_al_South_America_sites)
data(Venkataraman_et_al)
data(Winterhalder)
data(Yu_et_al)
data(Ziker)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\format{
    \enumerate{
        \item trip_id : Identifier for trip
        \item trip_id_orig : Identifier in raw data (for reference and cross-check)
        \item observed : 0/1 indicator of whether trip was directory observed or verbally reported
        \item trip_date : Calendar date on which trip occurred
        \item trip_year : Calendar year of trip
        \item julian_date_s : Standardized Julian date of trip
        \item day_id : Identifier for unique day in sample
        \item trip_duration : Duration in hours
        \item day_trip : 0/1 indicator of whether was a single day trip (1) or multiple days with camping (0)
        \item group_type : \code{Independent}, \code{Cooperative}, or \code{Pooled} foraging trips
        \item pooled : 0/1 indicator for pooled trips
        \item harvest : Harvest in kilograms
        \item forager_id : Unique indentifier for forager
        \item forager_id_orig : ID in raw data (for reference and cross-check)
        \item age_type : Either \code{Exact} age or \code{Interval} age
        \item age_dist_1 : When age_type \code{Exact}, age in years. Otherwise lower bound of interval.
        \item age_dist_2 : When age_type \code{Exact}, standard error of age. Otherwise upper bound of interval.
        \item dogs : 0/1 indicator of whether dogs were used on trip
        \item gun : 0/1 indicator of firearm use
        \item a_n_x : Variables referencing assistant number \code{n}, differing by availability in each data set
    }
}
\value{
}
\references{
Alvard: Alvard M. 1995. Shotguns and sustainable hunting in the Neotropics. Oryx 29:58-66.

Beckerman: Beckerman S, Lizarralde R. 2013. The Ecology of the Bar’: Rainforest Horticulturalists of South America. University of Texas Press.

Bird_Bird_Codding: Bird DW, Bird RB, Codding BF. 2009. In pursuit of mobile prey: Martu hunting strategies and archaeofaunal interpretation. American Antiquity 74:3-29.

Coad: Coad LM. 2008. Bushmeat hunting in Gabon: socio-economics and hunter behaviour. Doctoral dissertation, University of Cambridge.

Duda: Reyes-Garc’a V, GuЏze M, D’az-Reviriego I, Duda R, Fern‡ndez-Llamazares з, Gallois S, Napitupulu L, Orta-Mart’nez M, PyhЉlЉ A. 2016. The adaptive nature of culture: a cross-cultural analysis of the returns of local environmental knowledge in three indigenous societies. Current Anthropology 57:761-784.

Ellen: Ellen R. 1996. Individual strategy and cultural regulation in Nuaulu hunting. In: Ellen R, Fukui K. (eds) Redefining Nature: Ecology, Culture and Domestication, pp. 597-635. Berg, Oxford.

Fernandez_Llamazares: Reyes-Garc’a V, GuЏze M, D’az-Reviriego I, Duda R, Fern‡ndez-Llamazares з, Gallois S, Napitupulu L, Orta-Mart’nez M, PyhЉlЉ A. 2016. The adaptive nature of culture: a cross-cultural analysis of the returns of local environmental knowledge in three indigenous societies. Current Anthropology 57:761-784.

Franzen: Franzen M. 2006. Evaluating the sustainability of hunting: a comparison of harvest profiles across three Huaorani communities. Environmental Conservation 33:36-45.

Gallois: Reyes-Garc’a V, GuЏze M, D’az-Reviriego I, Duda R, Fern‡ndez-Llamazares з, Gallois S, Napitupulu L, Orta-Mart’nez M, PyhЉlЉ A. 2016. The adaptive nature of culture: a cross-cultural analysis of the returns of local environmental knowledge in three indigenous societies. Current Anthropology 57:761-784.

Gueze: Reyes-Garc’a V, GuЏze M, D’az-Reviriego I, Duda R, Fern‡ndez-Llamazares з, Gallois S, Napitupulu L, Orta-Mart’nez M, PyhЉlЉ A. 2016. The adaptive nature of culture: a cross-cultural analysis of the returns of local environmental knowledge in three indigenous societies. Current Anthropology 57:761-784.

Headland: Headland T. 1986. Why foragers do not become farmers: a historical study of a changing ecosystem and its effects on a Negrito hunter-gatherer group in the Philippines. Doctoral dissertation, University of Michigan.

Healey: Healey CJ. 1990. Maring Hunters and Traders: Production and Exchange in the Papua New Guinea Highlands. University of California Press.
 
Hill_Kintigh: Hill K, Kintigh K. 2009. Can anthropologists distinguish good and poor hunters? Implications for hunting hypotheses, sharing conventions, and cultural transmission. Current Anthropology 50:369-378.

Koster: Koster JM. 2008. Hunting with dogs in Nicaragua: An optimal foraging approach. Current Anthropology 49:935-944.

Kramer_Greaves: Kramer KL, Greaves RD. 2017. Why PumЋ Foragers retain a hunting and gathering way of life. In: Reyes-Garc’a V, PyhЉlЉ A (eds) Hunter-gatherers in a Changing World, pp. 109-126. Springer, Cham.

Lupo_Schmitt: Lupo KD, Schmitt DN. 2002. Upper Paleolithic net-hunting, small prey exploitation, and women's work effort: a view from the ethnographic and ethnoarchaeological record of the Congo Basin. Journal of Archaeological Method and Theory 9:147-179.

Marks: Marks SA. 1977. Hunting behavior and strategies of the Valley Bisa in Zambia. Human Ecology 5:1-36.

Napitupulu: Reyes-Garc’a V, GuЏze M, D’az-Reviriego I, Duda R, Fern‡ndez-Llamazares з, Gallois S, Napitupulu L, Orta-Mart’nez M, PyhЉlЉ A. 2016. The adaptive nature of culture: a cross-cultural analysis of the returns of local environmental knowledge in three indigenous societies. Current Anthropology 57:761-784.

Nielsen: Nielsen MR. 2006. Importance, cause and effect of bushmeat hunting in the Udzungwa Mountains, Tanzania: Implications for community based wildlife management. Biological Conservation 128:509-516.

Pacheco: Pacheco-Cobos L, Winterhalder B. 2015. Economic Benefits of Hunting with Dogs in the Context of Tropical Horticulture. 80th Annual Meeting of the Society for American Archaeology. San Francisco CA, USA

Pangau_Adam: Pangau-Adam M, Noske R, Muehlenberg M. 2012. Wildmeat or bushmeat? Subsistence hunting and commercial harvesting in Papua (West New Guinea), Indonesia. Human Ecology 40:611-621.
Ready   Ready E. 2016. Food, sharing, and social structure in an Arctic mixed economy. Doctoral dissertation, Stanford University.

Reyes-Garcia: Reyes-Garc’a V, GuЏze M, D’az-Reviriego I, Duda R, Fern‡ndez-Llamazares з, Gallois S, Napitupulu L, Orta-Mart’nez M, PyhЉlЉ A. 2016. The adaptive nature of culture: a cross-cultural analysis of the returns of local environmental knowledge in three indigenous societies. Current Anthropology 57:761-784.

Ross: Ross CT, Winterhalder B. 2018. Evidence for encounter-conditional, area-restricted search in a preliminary study of Colombian blowgun hunters. PLOS One, 13:e0207633.

Sillitoe: Sillitoe P. 2004. Managing Animals in New Guinea: Preying the Game in the Highlands. Routledge.

Siren: SirЋn AH, Wilkie DS. 2016. The effects of ammunition price on subsistence hunting in an Amazonian village. Oryx 50:47-55.

Trumble_Gurven: Gurven M, Kaplan H, Gutierrez M. 2006. How long does it take to become a proficient hunter? Implications for the evolution of extended development and long life span. Journal of Human Evolution 51:454-470.

Van_Vliet_et_al_Baego: van Vliet N, Vanegas L, Sandrin F, Cornelis D, Le Bel S, Dominique E, Gevais OO, Gaidet N, Fargeot C, Essiane E, Sicard J, Gely M, Lescuyer G, Billand A, Nasi R, Jepang C, Ayaya I, Broussolle L, Muhindo J, Houngbegnon F, Fagot R, Arnaud S, Zibe MB, Granier E, Kidiba A, Yapi F, Ngohouani D. 2017. Diagnostic approfondi pour la mise en oeuvre de la gestion communautaire

Van_Vliet_et_al_Djoutou: van Vliet N, Vanegas L, Sandrin F, Cornelis D, Le Bel S, Dominique E, Gevais OO, Gaidet N, Fargeot C, Essiane E, Sicard J, Gely M, Lescuyer G, Billand A, Nasi R, Jepang C, Ayaya I, Broussolle L, Muhindo J, Houngbegnon F, Fagot R, Arnaud S, Zibe MB, Granier E, Kidiba A, Yapi F, Ngohouani D. 2017. Diagnostic approfondi pour la mise en oeuvre de la gestion communautaire

Van_Vliet_et_al_Gabon: Van Vliet N, Nasi R. 2008. Hunting for livelihood in northeast Gabon: patterns, evolution, and sustainability. Ecology and Society 13.

Van_Vliet_et_al_Ingolo: van Vliet N, Vanegas L, Sandrin F, Cornelis D, Le Bel S, Dominique E, Gevais OO, Gaidet N, Fargeot C, Essiane E, Sicard J, Gely M, Lescuyer G, Billand A, Nasi R, Jepang C, Ayaya I, Broussolle L, Muhindo J, Houngbegnon F, Fagot R, Arnaud S, Zibe MB, Granier E, Kidiba A, Yapi F, Ngohouani D. 2017. Diagnostic approfondi pour la mise en oeuvre de la gestion communautaire

Van_Vliet_et_al_Ngombe: van Vliet N, Vanegas L, Sandrin F, Cornelis D, Le Bel S, Dominique E, Gevais OO, Gaidet N, Fargeot C, Essiane E, Sicard J, Gely M, Lescuyer G, Billand A, Nasi R, Jepang C, Ayaya I, Broussolle L, Muhindo J, Houngbegnon F, Fagot R, Arnaud S, Zibe MB, Granier E, Kidiba A, Yapi F, Ngohouani D. 2017. Diagnostic approfondi pour la mise en oeuvre de la gestion communautaire

Van_Vliet_et_al_Ovan: van Vliet N, Vanegas L, Sandrin F, Cornelis D, Le Bel S, Dominique E, Gevais OO, Gaidet N, Fargeot C, Essiane E, Sicard J, Gely M, Lescuyer G, Billand A, Nasi R, Jepang C, Ayaya I, Broussolle L, Muhindo J, Houngbegnon F, Fagot R, Arnaud S, Zibe MB, Granier E, Kidiba A, Yapi F, Ngohouani D. 2017. Diagnostic approfondi pour la mise en oeuvre de la gestion communautaire

Van_Vliet_et_al_Phalanga: van Vliet N, Vanegas L, Sandrin F, Cornelis D, Le Bel S, Dominique E, Gevais OO, Gaidet N, Fargeot C, Essiane E, Sicard J, Gely M, Lescuyer G, Billand A, Nasi R, Jepang C, Ayaya I, Broussolle L, Muhindo J, Houngbegnon F, Fagot R, Arnaud S, Zibe MB, Granier E, Kidiba A, Yapi F, Ngohouani D. 2017. Diagnostic approfondi pour la mise en oeuvre de la gestion communautaire

Van_Vliet_et_al_South_America_sites: van Vliet N, Mesa MP, Cruz-Antia D, de Aquino LJ, Moreno J, Nasi R. 2014. The uncovered volumes of bushmeat commercialized in the Amazonian trifrontier between Colombia, Peru & Brazil. Ethnobiology and Conservation 3.

Venkataraman_et_al: Endicott K. 1979. The hunting methods of the Batek negritos of Malaysia: A problem of alternatives. Canberra Anthropology 2:7-22.

Winterhalder: Winterhalder B. 1983. Boreal foraging strategies. In: Steegmann A.T. (eds) Boreal Forest Adaptations, pp. 201-241. Springer, Boston, MA

Yu_et_al: Ohl-Schacherer J, Shepard Jr GH, Kaplan H, Peres CA, Levi T, Yu DW. 2007. The sustainability of subsistence hunting by Matsigenka native communities in Manu National Park, Peru. Conservation Biology 21:1174-1185.

Ziker: Ziker JP. 2002. Peoples of the tundra: Northern Siberians in the post-communist transition. Waveland Press.
}
\seealso{}
\examples{
data(Alvard)
head(Alvard)
plot( harvest ~ trip_duration , Alvard )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }


# ecoToxic.R 
# Part of the SEHmodel package.
#
# Copyright (C) 2015        Melen Leclerc <melen.leclerc@rennes.inra.fr>
#                           Jean-Francois Rey <jean-francois.rey@paca.inra.fr>
#                           Samuel Soubeyrand <Samuel.Soubeyrand@avignon.inra.fr>
#                           Emily Walker <emily.walker@avignon.inra.fr>
#                           INRA - BioSP Site Agroparc - 84914 Avignon Cedex 9
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

#' @title Method ecoToxicological
#' @description Generic method on Landscape and Individuals objects for ecotoxicological.
#' @name ecoToxic
# @param objectL A Landscape object
# @param objectI An Individuals object
#' @param ... other parameters
#' @rdname ecoToxic
#' @exportMethod ecoToxic
setGeneric(name="ecoToxic",
           def=function(objectL,objectI,...)
             standardGeneric("ecoToxic")
)

#' ecoToxic
#' @description Internal concentration of contaminants within individuals
# @details ecotoxic details
#' @param objectL A Landscape object
#' @param objectI An Individuals object
#' @param objectT A 3D array of Toxic Dispersion over time [t,x,y], first indice is time 
#' @param mintime Time to start simulation (default = 1)
#' @param maxtime Time to end simulation
#' @param kin \% of toxic ingestion
#' @param kout \% of toxic elimination
#' @param deltat \% of a time step
#' @return An Individuals object
#' @aliases ecotoxicological-method
#' @rdname ecoToxic
setMethod(f="ecoToxic",
          signature=c("Landscape","Individuals"),
          definition=function(objectL,objectI,objectT,mintime=1,maxtime,kin=0.25,kout=0.5,deltat=0.1) {
            xstep=(objectL@xmax-objectL@xmin)/length(objectT[mintime,,1])
            ystep=(objectL@ymax-objectL@ymin)/length(objectT[mintime,1,])
            for(i in 1:objectI@n) {
              #cmilieu=array(NA,length(objectT[,1,1]))
              cmilieu<-objectT[,round((objectI@coordinate[i]@coords[1]-objectI@xmin)/xstep),round((objectI@coordinate[i]@coords[2]-objectI@ymin)/ystep)]
              if(mintime > objectI@dob[i]) { cint<-objectI@intern_toxic[i,mintime-1] }
              else {cint<-0}
              conctemp<-conc.int(cmil=cmilieu,cint_start=cint,kin=kin,kout=kout,min.time=max(mintime,objectI@dob[i]),min(maxtime,objectI@dob[i]+objectI@life_duration[i]-1),deltat=deltat)
              objectI@intern_toxic[i,max(mintime,objectI@dob[i]):min(maxtime,(objectI@dob[i]+objectI@life_duration[i]-1))]<-conctemp
            }
            return(objectI)
        }
)



#' @title Plot internal toxic concentration
#'
#' @description Plot a time serie of internal toxic concentration for a given individual.
#' @name plotEcotoxic
#' @param objectL A Landscape object
#' @param objectI An Individuals object
#' @param objectT A 3d array of Toxic intensity over the time [t,x,y], (first indice is time) see \code{\link{toxicIntensity}}
#' @param numind An individual id
# @rdname ecoToxic
#' @export
plotEcoToxic<-function(objectL,objectI,objectT,numind=8) {
  mintime=1
  xstep=(objectL@xmax-objectL@xmin)/length(objectT[mintime,,1])
  ystep=(objectL@ymax-objectL@ymin)/length(objectT[mintime,1,])
  par(mar=c(5.1,4.1,4.1,2.1))
  for(i in numind ) {
    plot(objectT[,objectI@coordinate[i]@coords[1]/xstep,objectI@coordinate[i]@coords[2]/ystep],col=2,xlab="time",ylab="toxic",type="b",pch=16)
    lines(objectI@intern_toxic[i,],col=1,type="b",pch=16)
    points(objectI@dob[i],0,col=3,pch=16,cex=1.5)
    legend("topright",c("internal concentration","local conc","date of birth"),col=c(1,2,3),pch=16 )
    title(main=paste("Individual #",numind))
  }
}

# @export
#getIndividualsInfo <- function(objectL,objectI,objectT) {
#  mintime=1
#  xstep=(objectL@xmax-objectL@xmin)/length(objectT[mintime,,1])
#  ystep=(objectL@ymax-objectL@ymin)/length(objectT[mintime,1,])
#  dfres<-data.frame()
#  for(i in 1:objectI@n) {
#    dfind<-objectI[i]
#    df<-cbind(dfind,"local_intensity"=matrix(objectT[,round((objectI@coordinate[i]@coords[1]-objectI@xmin)/xstep),round((objectI@coordinate[i]@coords[2]-objectI@ymin)/ystep)],ncol=length(objectT[,round((objectI@coordinate[i]@coords[1]-objectI@xmin)/xstep),round((objectI@coordinate[i]@coords[2]-objectI@ymin)/ystep)])))
#    dfres<-rbind(dfres,df)
#  }
#  return(dfres)
#}

###############
### PRIVATE ###
###############

# emily : modele ecotoxico simple : ashauer et al 2007
# cint : concentration interne à suivre dans le temps + seuil de mortalité à fixer sur cette concentration
# cmil : concentration dans le milieu (vecteur de 30 valeurs pour 30 jours)
# cint_start : concentration interne au debut
# kin : constante d'absorption ( 25% du pollen d'un pixel est absorbé par jour par larve)
# kout : constante d'élimination (50% du pollen est éliminé d'un jour au suivant)
# min.time temps de début
# max.time temps de fin
# deltat 
conc.int=function(cmil,cint_start=0,kin=0.25,kout=0.5,min.time=1,max.time,deltat=0.01) {
  if(min.time > max.time) {stop("ERROR : ecoToxic time value error")}
  cint=numeric(0)  
  cint_previous=cint_start
  for (t in min.time:max.time) {
    cmilieut=cmil[t]
    cint_temp=cint_previous
    j=1
    for (dt in seq(0,0.99,deltat)) { # max.time+1
      # attention vecteur indices=entiers
      cint_temp[j+1]<-(cint_temp[j] + deltat*kin*cmilieut) / (1 + deltat*kout)
      j=j+1
    }
    cint[t]<-cint_temp[length(cint_temp)]
    cint_previous<-cint[t]
  }
  return(cint[min.time:max.time])
}

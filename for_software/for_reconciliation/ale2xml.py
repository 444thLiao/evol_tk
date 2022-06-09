"""
Script modified from the https://github.com/WandrilleD/recPhyloXML/tree/master/python3

By thliao 20220531
"""

import ete3
from ete3 import Tree


# default definitions/tags of recphyloXML
RECPHYLOTAG = "recPhylo"
RECTREETAG = "recGeneTree"
SPTREETAG = "spTree"

EVENTTAGCORRESPONDANCE = {
    "D": "duplication",
    "S": "speciation",
    "C": "leaf",
    "L": "loss",
    "Bo": "bifurcationOut",
    "bro": "branchingOut",
    "Tb": "transferBack",
    "SL": "speciationLoss",
    "broL": "branchingOutLoss",
}

## helper function to get XML lines for a simple XML tree. typically used for the species tree here
def myBasicTreeXMLLinesAux(tree):
    """
    Takes:
        - tree (ete3.TreeNode)
    Returns:
        (list): list of xml lines With names of internal node 
    """

    indentChar = "  "
    lines = ["<clade>"]
    lines.append(indentChar + "<name>" + tree.name + "</name>")
    for c in tree.children:
        tmp = myBasicTreeXMLLinesAux(c)
        for l in tmp:
            lines.append(indentChar + l)
    lines.append("</clade>")
    return lines

def myBasicTreeXMLLines(tree):
    """
    Takes:
        - tree (ete3.TreeNode)
    Returns:
        (list): list of xml lines
    """
    lines = ["<phylogeny>"]
    indentChar = "  "
    tmp = myBasicTreeXMLLinesAux(tree)
    for l in tmp:
        lines.append(indentChar + l)
    lines.append("</phylogeny>")
    return lines

class RecEvent:
    def __init__(self, eventCode, species, ts=None, additionnalInfo={}):
        """
        Takes:
            - eventCode (str) : a code indicating the recEvent event
            - species (~) : a identifier for the specie the event takes place in
            - ts (int or None) [default= None] : the time slice the events happens at, if applicable
            - additionnalInfo (dict) [default= {}] : keys are expected to be some property tag and values the associated information
        """
        self.eventCode = eventCode
        self.species = species
        self.timeSlice = ts
        self.additionnalInfo = additionnalInfo.copy()

    def __str__(self):
        eventName = self.eventCode
        #        if EVENTTAGCORRESPONDANCE.has_key(self.eventCode):
        if self.eventCode in EVENTTAGCORRESPONDANCE:
            eventName = EVENTTAGCORRESPONDANCE[self.eventCode]
        L = [str(self.eventCode), "spe=" + str(self.species)]
        if not self.timeSlice is None:
            L.append("ts=" + str(self.timeSlice))
        for k, v in self.additionnalInfo.items():
            L.append(str(k) + "=" + str(v))
        return " ".join(L)
    def nwkstr(self):
        """ tmp simplistic version """
        s = str(self.species)
        s += "."

        s += str(self.eventCode)
        return s
    def makeRecXMLstr(self, speciesNames):
        if self.eventCode == "N":
            return ""
        eventName = self.eventCode
        #        if EVENTTAGCORRESPONDANCE.has_key(self.eventCode):
        if self.eventCode in EVENTTAGCORRESPONDANCE:
            eventName = EVENTTAGCORRESPONDANCE[self.eventCode]
        spe = str(self.species)
        if spe in speciesNames.keys():
            spe = speciesNames[spe]
        S = "<" + str(eventName)
        if self.eventCode != "Bo":
            S += " "
            if self.eventCode == "Tb":
                S += "destinationSpecies="
            else:
                S += "speciesLocation="
            S += '"' + str(spe) + '"'

            if not self.timeSlice is None:
                S += " ts=" + '"' + str(self.timeSlice) + '"'
        propertyName = "confidence"
        #        if self.additionnalInfo.has_key(propertyName):
        if propertyName in self.additionnalInfo:
            S += (
                " "
                + propertyName
                + "="
                + '"'
                + self.additionnalInfo[propertyName]
                + '"'
            )
        if self.eventCode == "C":
            propertyName = "geneName"
            #            if self.additionnalInfo.has_key(propertyName):
            if propertyName in self.additionnalInfo:
                S += (
                    " "
                    + propertyName
                    + "="
                    + '"'
                    + self.additionnalInfo[propertyName]
                    + '"'
                )
        S += ">"
        S += "</" + eventName + ">"
        # print self.eventCode , "->", S
        return S

class ReconciledTree(ete3.TreeNode):
    def __init__(self):
        ete3.TreeNode.__init__(self)

        self.name = ""
        self.eventRecs = []

    def setName(self, name):
        """ NB : name can be any object with a __str__() fc (like an int for instance) """
        self.name = name

    def getEvents(self):
        return self.eventRecs

    def getEvent(self, i):
        return self.eventRecs[i]

    def popEvent(self, i):
        return self.eventRecs.pop(i)

    def addEvent(self, e, append=True):
        """
            Takes:
                - e (RecEvent) : the reconciliation event to add
                - append (bool) [default=True] : if True adds the event at the end of the list else, adds it at the beginning
        """
        if append:
            self.eventRecs.append(e)
        else:
            self.eventRecs.insert(0, e)
        return

    def getTreeStrAux(self):

        L = []
        L.append("name : " + str(self.name))
        L.append("events :")
        for e in self.eventRecs:
            L.append("  " + str(e))
        ChL = []
        for c in self.get_children():
            ChL += ["  " + s for s in c.getTreeStrAux()]
        if len(ChL) > 0:
            L.append("children :")
            L += ChL
        return L

    def getTreeStr(self):
        return "\n".join(self.getTreeStrAux())

    def getTreeNewickAux(self, sep="|", topoOnly=False):
        s = ""

        ChL = []
        for c in self.get_children():
            ChL.append(c.getTreeNewickAux(sep, topoOnly))
        if len(ChL) > 0:
            s += "("
            s += ",".join(ChL)
            s += ")"

        s += str(self.name)
        if not topoOnly:
            s += sep
            s += sep.join([e.nwkstr() for e in self.eventRecs])

        return s

    def getTreeNewick(self, sep="|", topoOnly=False):
        return self.getTreeNewickAux(sep, topoOnly) + ";"

    def getTreeRecPhyloXMLAux(self, speciesNames={}, topoOnly=False):

        L = []
        L.append("<clade>")
        L.append("  <name>" + str(self.name) + "</name>")
        if not topoOnly:
            L.append("  <eventsRec>")
            for e in self.eventRecs:
                s = e.makeRecXMLstr(speciesNames)
                if s == "":
                    continue
                L.append("    " + s)
            L.append("  </eventsRec>")
        ChL = []
        for c in self.get_children():
            ChL += ["  " + s for s in c.getTreeRecPhyloXMLAux(speciesNames, topoOnly)]
        if len(ChL) > 0:
            L += ChL
        L.append("</clade>")
        return L

    def getTreeRecPhyloXML(self, speciesNames={}, topoOnly=False):
        Lines = self.getTreeRecPhyloXMLLines(speciesNames, topoOnly)
        return "\n".join(Lines)

    def getTreeRecPhyloXMLLines(self, speciesNames={}, topoOnly=False):
        Lines = ["<recGeneTree>"]
        Lines.append('<phylogeny rooted="true">')
        tmp = self.getTreeRecPhyloXMLAux(speciesNames, topoOnly)
        for l in tmp:
            Lines.append("    " + l)
        Lines.append("</phylogeny>")
        Lines.append("</recGeneTree>")
        return Lines

    def countEvents(self):
        """

        Returns:
            (dict) : keys are recPhyloXML event tags, values are the number of times these events occur in the tree
        """
        devent = {}
        for e in self.eventRecs:
            code = e.eventCode
            #            if not devent.has_key(code):
            if not code in devent:
                devent[code] = 0
            devent[code] += 1

        for c in self.get_children():
            tmp = c.countEvents()
            for k in tmp.keys():
                #                if not devent.has_key(k):
                if not k in devent:
                    devent[k] = 0
                devent[k] += tmp[k]

        return devent

    def getEventsSummary(
        self,
        speciesTree,
        includeTransferReception=True,
        includeTransferDeparture=False,
        speciesIdFeature="name",
    ):
        """
        *recursive function*

        Takes:
             - speciesTree (ete3.Tree) : the species tree used for the reconciliation, necessary to assign a species to loss events.
             - includeTransferReception [default = True]  : Whether or not to includes events of  TransferReception (transferBack tag) in the counts.
             - includeTransferDeparture [default = False] : Whether or not to includes events of  TransferDeparture (branchingOut tag) in the counts.
             - speciesIdFeature (str) [default = "name"] : the feature to use as Id in the species tree (by default the name is used)

        Returns:
            (dict):
                    keys are events type among : "duplication" , "loss" , "transferReception" , "transferDeparture"
                    values are  lists of species id

        """

        EventsSummary = {"duplication": [], "loss": []}

        if includeTransferReception:
            EventsSummary["transferReception"] = []

        if includeTransferDeparture:
            EventsSummary["transferDeparture"] = []

        for i, e in enumerate(self.eventRecs):

            evtCode = e.eventCode

            species = e.species

            report = False

            reportTD = False
            reportTR = False

            if (
                evtCode
                in (EVENTTAGCORRESPONDANCE["SL"], EVENTTAGCORRESPONDANCE["broL"])
            ) or (evtCode in ("SL", "broL")):
                species = self.getLostSpecies(i, speciesTree, speciesIdFeature)
                evtCode = "loss"
                report = True

            elif evtCode == EVENTTAGCORRESPONDANCE["D"] or evtCode == "D":
                evtCode = "duplication"
                report = True

            elif evtCode == EVENTTAGCORRESPONDANCE["L"] or evtCode == "L":
                evtCode = "loss"
                report = True

            if report:  ##reporting duplication or loss
                EventsSummary[evtCode].append(species)

            if includeTransferReception and (
                evtCode == EVENTTAGCORRESPONDANCE["Tb"] or evtCode == "Tb"
            ):
                evtCode = "transferReception"
                EventsSummary[evtCode].append(species)

            elif includeTransferDeparture and (
                (
                    evtCode
                    in [EVENTTAGCORRESPONDANCE["bro"], EVENTTAGCORRESPONDANCE["broL"]]
                )
                or (evtCode in ["bro", "broL"])
            ):
                evtCode = "transferDeparture"
                EventsSummary[evtCode].append(species)

        for c in self.get_children():

            tmp = c.getEventsSummary(
                speciesTree,
                includeTransferReception,
                includeTransferDeparture,
                speciesIdFeature,
            )

            for k, v in tmp.items():
                EventsSummary[k] += v

        return EventsSummary

    def sameSpeciesAsParent(self, parent=None):
        """ returns True if the first event of the node has the same species as the last event of its parent , False otherwise (and if self is the root)

            if the parent is given, is it used, otherwise we look for it in the structure
        """

        if parent is None:
            if self.is_root():
                return False
            parent = self.up

        lastParentSp = parent.getEvents()[-1].species

        firstSp = self.getEvents()[0].species

        return firstSp == lastParentSp

    def getLostSpecies(self, evtIndex, speciesTree, speciesIdFeature="name"):
        """
        given the index of an event of *Loss (speciationLoss for instance) in this nodes,
        this function returns the id of the species where the loss occured
        (this function is useful because speciationLoss event references the species of the speciation rather than the species of the loss)

        Takes:
            - evtIndex (int) : index of the loss event whose lost species we want to know
            - speciesTree (ete3.Tree) : the species tree used for the reconciliation
            - speciesIdFeature (str) [default = "name"] : the feature to use as Id in the species tree (by default the name is used)

        Returns:
            (str) : id of the species where the loss occured
            or
            None : if the indicated event does not correspond to a loss
        """

        evtCode = self.getEvent(evtIndex).eventCode

        if not evtCode in [
            "SL",
            "broL",
            EVENTTAGCORRESPONDANCE["SL"],
            EVENTTAGCORRESPONDANCE["broL"],
        ]:
            return None

        if evtCode in ("broL", EVENTTAGCORRESPONDANCE["broL"]):
            return self.getEvent(
                evtIndex
            ).species  ## in branchingOutLoss, the species of the lost lineage is the same as the one of the branchingOut

        ## We know this is a speciationLoss event.
        ## Note that speciationLoss events are NEVER the last event in eventsRec (as they are neither a bifurcation nor a leaf event).
        ## so we should be able to safely ask for the event after the current one.

        lostSpeciesSister = self.getEvent(evtIndex + 1).species

        Speciesnode = speciesTree.search_nodes(**{speciesIdFeature: lostSpeciesSister})

        if len(Speciesnode) != 1:
            raise Error(
                "error:",
                len(Speciesnode),
                "with Id",
                lostSpeciesSister,
                "(1 expected).",
            )

        lostSpeciesNode = Speciesnode[0].get_sisters()[0]

        lostSpeciesId = getattr(lostSpeciesNode, speciesIdFeature)

        return lostSpeciesId

class ReconciledTreeList:
    """
    This object represents a group of reconciled tree.
    Usually they would be reconciled with the same species tree, which can also be added to this object.

    Atributes:
        - self.spTree   : the species tree these trees are reconciled with (or None if no species tree is specified)
        - self.recTrees : the reconciled trees

    """

    def __init__(self, spTree=None, recTrees=[]):
        """
        Takes:
            spTree (ete3.Tree) [ default = None ] : facultative species tree
            recTrees (list) [ default = [] ] : list of ReconciledTree instance (see above)
        """

        self.spTree = spTree
        self.recTrees = recTrees[:]

    def setSpTree(self, ST):
        """
        Simply sets a trees as the object species tree.

        Takes:
            - ST (ete3.Tree) : a species tree
        """
        self.spTree = ST

    def append(self, RT):
        """
        Appends a reconciled tree to the object.

        Takes:
            - RT (ReconciledTree) : a reconciled tree
        """
        self.recTrees.append(RT)

    def __getitem__(self, i):
        """
        returns a reconciled tree from the object.

        Takes:
            - i (int) : index of the desired reconciled tree

        Returns:
            (ReconciledTree) : the reconciled tree at the desired index
            OR IndexError if the index is invalid (ie. too high)
        """

        if i >= len(self.recTrees):
            raise IndexError(
                "Index out of range. There are no reconciled tree with index "
                + str(i)
                + "."
            )

        return self.recTrees[i]

    def __len__(self):
        """
        Returns:
            (int) : number of ReconciledTree in this instance
        """

        return len(self.recTrees)

    def hasSpTree(self):
        """
        Returns:
            (bool) : True if there this instance has a species tree (ie. self.spTree is not None), False otherwise
        """
        return not self.spTree is None

    def getRecPhyloXMLLines(self):
        """
        Returns:
            (list) : list of lines of the recPhyloXML representation of this object
                     (NB : the lines do not have a '\n' at their end.)
        """
        lines = []
        lines.append("<" + RECPHYLOTAG + ">")

        offset = 1
        offsetChar = "  "

        if self.hasSpTree():
            lines.append(offsetChar * offset + "<" + SPTREETAG + ">")
            offset += 1
            spLines = myBasicTreeXMLLines(self.spTree)
            for l in spLines:
                lines.append(offsetChar * offset + l)
            offset -= 1
            lines.append(offsetChar * offset + "</" + SPTREETAG + ">")

        for RT in self.recTrees:
            recLines = RT.getTreeRecPhyloXMLLines()
            for l in recLines:
                lines.append(offsetChar * offset + l)

        lines.append("</" + RECPHYLOTAG + ">")

        return lines

    def getEventsSummary(
        self,
        includeTransferReception=True,
        includeTransferDeparture=False,
        speciesIdFeature="name",
        indexBySpecies=False,
    ):
        """
        Retrieve an event summary over all the trees in the object
        !!only works if there is a species tree assigned to the object!!


        Takes:
             - includeTransferReception [default = True]  : whether or not to includes events of  TransferReception (transferBack tag) in the counts.
             - includeTransferDeparture [default = False] : whether or not to includes events of  TransferDeparture (branchingOut tag) in the counts.
             - speciesIdFeature (str) [default = "name"] : the feature to use as Id in the species tree (by default the name is used)
             - indexBySpecies (str) [default = False] : if True, the returned dictionnary will have species as keys and event counts as values.

        Returns:
            (dict):
                    keys are events type among : "duplication" , "loss" , "transferReception" , "transferDeparture"
                    values are  lists of species id

                   OR, if indexBySpecies=True:
                       keys are species id
                       values are dict with keys among "duplication" , "loss" , "transferReception" , "transferDeparture"
                                            and values as counts of the events in each species

        """

        if not self.hasSpTree():
            raise Error(
                "error : can't get an events summary when no species tree has been assigned."
            )

        EventsSummary = {"duplication": [], "loss": []}

        if includeTransferReception:
            EventsSummary["transferReception"] = []

        if includeTransferDeparture:
            EventsSummary["transferDeparture"] = []

        for RT in self.recTrees:
            tmp = RT.getEventsSummary(
                self.spTree,
                includeTransferReception,
                includeTransferDeparture,
                speciesIdFeature,
            )

            for k, v in tmp.items():
                EventsSummary[k] += v

        if indexBySpecies:

            tmp = {}

            for n in self.spTree.traverse():

                tmp[getattr(n, speciesIdFeature)] = {}

            for e in EventsSummary.keys():

                for sp in tmp.keys():
                    tmp[sp][e] = 0

                for sp in EventsSummary[e]:
                    tmp[sp][e] += 1

            EventsSummary = tmp

        return EventsSummary

def completeTreeNames(tree, useBS=False):
    """
    Takes:
        - tree (ete3.Tree)
        - useBS (bool) [default = False] : uses bootstrap to name nodes

    Returns:
        (ete3.Tree) : the tree, but where the nodes without a name now have one that correspond
                    to their post-order OR their bootstrap

    """

    for i, n in enumerate(tree.traverse("postorder")):
        if n.name == "":
            print(n.support)
            if useBS:
                n.name = str(int(n.support))
            else:
                n.name = str(i)
    return tree


# def treatLossEvent(node, LossEvent , keptChildNameSuffix = ".c"):
#    """
#    Takes:
#        - node (ReconciledTree)
#        - LossEvent (recEvent)
#        - keptChildNameSuffix (str) [default = ".c"] : suffix to add to the name of the new child of node that is NOT a loss
#    """
#
#
#    # 1. create the loss child
#
#    #determining ts if there is one
#    lostTS = LossEvent.timeSlice
#    if not lostTS is None:
#        lostTS -= 1 # one TS more recent than the parent
#
#    lossNode = ReconciledTree()
#    lossNode.addEvent( RecEvent("loss" , "", ts= lostTS ) )
#    lossNode.name="LOSS"
#
#    # 2. create the kept child
#
#    keptNode = ReconciledTree()
#    keptNode.name = node.name + keptChildNameSuffix
#
#    # 4. branching loss and kept to original node
#
#    node.add_child(lossNode)
#    node.add_child(keptNode)
#
#    # 5. editing the event and adding to node
#
#    e = LossEvent.eventCode
#    LossEvent.eventCode = e.rpartition("L")[0]
#
#    node.addEvent(LossEvent)
#
#    return keptNode


def parse_node_annotation(node_annotation, isLeaf=False, isDead=False, isUndated=False):
    """
    Takes:
        - node_annotation (str): reconciliation annotation on a node

    Returns:
        (list): list of dicts coding a particular event
    """
    # print("annot : ",node_annotation , isUndated)

    l_events = []

    if len(node_annotation) != 0:

        if node_annotation.startswith("."):
            node_annotation = node_annotation[1:]
        tmp_ann = node_annotation.split(".")

        ##further splitting multiple transfer
        s_ann = []
        for ann in tmp_ann:
            if ann.count("@") < 1:
                s_ann.append(ann)
                continue
            ## in case of transfer and loss like: T@27|SYNP2@26|ACAM1
            new_anns = ann.split("@")

            s_ann.append("@".join(new_anns[0:2]))  ##first tranfer, a transfer out

            for a in new_anns[2:]:  ##for each transfer after that (should be only one)
                s_ann.append("@" + a)

        for ann in s_ann:
            if len(ann) == 0:
                raise Exception("empty annotation")

            if ann[0].isdigit():  ##starts with a number spe,dup or spe+loss
                if ann.isdigit():  ##only numbers: spe or spe+loss
                    target = ann

                    ts = int(target)

                    if isUndated:
                        ts = None

                    l_events.append(RecEvent("S", target, ts))
                    continue

            if ann.startswith("T@"):  ##Transfer out

                ## time slice of the source
                source_ts = None
                source_sp = None

                if isUndated:
                    ## of the shape : "T@D->A"
                    source_sp = ann[2:].partition("->")[0]
                else:
                    source_ts, junk, source_sp = ann[2:].partition(
                        "|"
                    )  ## partitionning something like T@22|22
                    source_ts = int(source_ts)

                ##adding the event
                l_events.append(RecEvent("bro", source_sp, source_ts))

            if ann.startswith("@"):  # or ann.startswith("Tb@"):##transfer in or back

                pre = 3  ##cropping size
                if ann.startswith("@"):
                    pre = 1

                target_ts, junk, target_sp = ann[pre:].partition(
                    "|"
                )  ## partitionning something like @22|22 to get the time slice and specie

                ##adding the event
                l_events.append(RecEvent("Tb", target_sp, int(target_ts)))

            if ann.startswith("Tb@"):
                l_events.append(RecEvent("Bo", "-1"))

            if ann.startswith("D@"):  ##Duplication

                ts = None
                sp = None
                if isUndated:
                    sp = ann[2:]
                else:
                    ts, junk, sp = ann[2:].partition(
                        "|"
                    )  ## partitionning something like D@23|22
                    ts = int(ts)

                l_events.append(RecEvent("D", sp, ts))

    if isLeaf and (len(l_events) == 0 or l_events[-1].eventCode != "C"):
        ts = 0
        if isUndated:
            ts = None
        l_events.append(
            RecEvent("C", None, ts)
        )  ##temporary placeholder for the leaf species

    if isDead:  ## we start in the dead so the first event MUST be Bout or Tloss

        if not l_events[0].eventCode in ["Tb", "Bo"]:

            target_ts = l_events[0].timeSlice
            target_sp = l_events[0].species

            l_events.insert(0, RecEvent("Tb", target_sp, target_ts))

    ##adding loss labels
    for i in range(len(l_events) - 1):  ##all events but the last one
        if l_events[i].eventCode in ["bro", "S"]:
            l_events[i].eventCode += "L"

    return l_events


def separateLeafNameFromLeafAnnotation(leafName, sepSp="_", sepAnnot=(".", "@")):
    """
    Takes:
        - leafName (str) : name of a leaf, potentially containing reconciliation information (exemple: "g_g3.T@4|3@1|g" )
        - sepSp (str) [default = "_" ] : separator between species name and gene name
        - sepAnnot (tuple) [default = (".","@") ] : possible separators between gene name and annotations

    Returns:
        (tuple)
            (str) : gene name
            (str) : reconciliation annotation  (empty string if there is none)

    """
    spName, j, gNameAndAnnot = leafName.partition(sepSp)

    x = 0
    AnnotFound = False

    while (not AnnotFound) and (x < len(gNameAndAnnot)):

        if gNameAndAnnot[x] in sepAnnot:
            AnnotFound = True
            break
        x += 1
    # print "->", leafName[:x] , leafName[x:]
    return spName + sepSp + gNameAndAnnot[:x], gNameAndAnnot[x:]


def getLeafSpeciesFromLeafName(leafName, sepSp="_"):
    """
    Takes:
         - leafName (str) : name of a leaf, in the format: species separator gene
         - sepSp (str) [default = "_" ] : separator between species name and gene name

    Returns:
        (str) : species name
    """

    return leafName.partition(sepSp)[0]


def ALEtreeToReconciledTree(ALEtree, isDead=False, isUndated=False, sepSp="_"):
    """
    Recursively builds the reconciled tree

    Takes:
        - ALEtree (ete3.TreeNode) : a tree read from a reconciled ttree in the ALE format (ie. reconciliation annotations are in the .name field)
        - isDead (bool) [default = False] : indicates whether or not this lineage starts in a dead/unsampled species

    Returns:
        (ReconciledTree)
    """

    isLeaf = ALEtree.is_leaf()

    annotation = None
    name = ALEtree.name
    if isLeaf:
        name, annotation = separateLeafNameFromLeafAnnotation(ALEtree.name, sepSp=sepSp)
        # print("leaf parsing :", name , annotation)
    else:
        annotation = ALEtree.name

    # print "name : ",ALEtree.name

    events = parse_node_annotation(
        annotation, isLeaf, isDead=isDead, isUndated=isUndated
    )

    if isLeaf:
        ## we specify the species of the leaf event
        events[-1].species = getLeafSpeciesFromLeafName(ALEtree.name, sepSp=sepSp)

    # print [str(e) for e in events]

    if events[-1].eventCode == "Bo":
        isDead = True  ## means that both children must begin by an event in the dead

    RT = ReconciledTree()
    RT.setName(name)

    current = RT

    for e in events:
        #        if e.eventCode.endswith("L"):
        #            #print "plep"
        #            current = treatLossEvent(current, e , ".c")
        #        else:
        #            current.addEvent(e)
        #
        current.addEvent(e)

    for c in ALEtree.children:  ##recursion on successors
        current.add_child(ALEtreeToReconciledTree(c, isDead, isUndated, sepSp=sepSp))

    return RT


def refineReconciledTreeWithTransferBack(RT):
    """
    adds transferBack events where they were omitted

    Takes:
        - RT (ReconciledTree) : a reconciled tree obtained from an ale string

    """

    for n in RT.traverse():
        if n.is_root():
            continue
        lastEventParent = n.up.getEvent(-1)
        if lastEventParent.eventCode in ["branchingOut", "bro"]:
            ## parent is an outgoing transfer
            firstEvent = n.getEvent(0)
            if firstEvent.species == lastEventParent.species:
                continue  ## this is the "kept" child --> continue
            if firstEvent.eventCode in ["Tb", "transferBack", "Bo", "bifurcationOut"]:
                continue  ## already the correct annotation
            TbEvent = RecEvent("Tb", species=firstEvent.species)
            n.addEvent(TbEvent, append=False)


def refineReconciledTreeLosses(RT, spTree):
    """
    adds species to the losses

    Takes:
        - RT (ReconciledTree) : a reconciled tree obtained from an ale string
        - spTree (ete3.Tree) : a species tree

    """

    for n in RT.traverse():
        if n.is_root():
            continue

        firstEvent = n.getEvent(0)

        if first.eventCode in ["L", "loss"]:
            ## loss!

            lastEventParent = n.up.getEvent(-1)

            if firstEvent.species == lastEventParent.species:
                continue  ## this is the "kept" child --> continue

            if firstEvent.eventCode in ["Tb", "transferBack", "Bo", "bifurcationOut"]:
                continue  ## already the correct annotation

            TbEvent = RecEvent("Tb", species=firstEvent.species)

            n.addEvent(TbEvent, append=False)


def ConvertRTtoLossIndepVersion(RT, speciesTree=None, keptChildNameSuffix=".c"):
    """
    *modifies RT*
    *RECURSIVE*

    Takes:
        - RT (ReconciledTree): reconciled tree or subtree to convert
        - speciesTree (ete3.Tree) [default = None] : species tree
        - keptChildNameSuffix (str) [default = ".c"] : suffix to add to the name of the new child of node that is NOT a loss
    """

    for i, e in enumerate(RT.eventRecs):

        if len(e.eventCode) > 1 and e.eventCode.endswith("L"):

            species = ""

            if not speciesTree is None:
                species = RT.getLostSpecies(i, speciesTree)

            lostTS = e.timeSlice
            if not lostTS is None:
                lostTS -= 1

            MakeLossIndependentNode(
                RT,
                i,
                lostSpecies=species,
                lostTS=lostTS,
                keptChildNameSuffix=keptChildNameSuffix,
            )

    for c in RT.children:
        ConvertRTtoLossIndepVersion(c, speciesTree, keptChildNameSuffix)

    return


def MakeLossIndependentNode(
    node,
    LossIndex,
    lostSpecies="",
    lostTS=None,
    lostAdditional={},
    keptChildNameSuffix=".c",
):
    """
    *modifies node*

    Takes:
         - node (ReconciledTree): reconciled node where the *Loss event occurs
         - LossIndex (int): index of the speciationLoss or branchingOutLoss event
         - lostSpecies (str) [default = ""] : species of the loss
         - lostTS (int) [default = None]: timeSlice is the loss
         - lostAdditional [default = {}]: additional information to give to the new loss event
         - keptChildNameSuffix (str) [default = ".c"] : suffix to add to the name of the new child of node that is NOT a loss
    """

    # print( MakeLossIndependentNode , node , lostTS )

    # 1. create the loss child

    lossNode = ReconciledTree()
    lossNode.addEvent(
        RecEvent("loss", lostSpecies, ts=lostTS, additionnalInfo=lostAdditional)
    )
    lossNode.name = "LOSS"

    # 2. create the kept child

    keptNode = ReconciledTree()
    keptNode.name = node.name + keptChildNameSuffix

    while len(node.eventRecs) > (LossIndex + 1):
        # print LossIndex, LossIndex+1 , len(node.eventRecs)
        keptNode.addEvent(node.popEvent(LossIndex + 1))

    # 3. link children to kept child

    while len(node.children) > 0:
        c = node.children[0]
        c.detach()
        keptNode.add_child(c)

    # 4. branching loss and kept to original node

    node.add_child(lossNode)
    node.add_child(keptNode)

    # 5. editing the event

    e = node.eventRecs[LossIndex].eventCode
    node.eventRecs[LossIndex].eventCode = e[:-1]

    return

TESTING_SPECIES_TREE_NEWICK = "(((a,b)2,(c,d)1)4,((e,f)3,g)5)6;"
TESTING_GENE_TREE_NEWICK = "((a_a1,b_b1@0|b).4.2.T@1|a,((((d_d1@0|d,e_e2.3)Tb@3|3,g_g2)T@4|g,g_g3.T@4|3@1|g).5,(e_e1.3,g_g1).5)D@5|5).6;"

# undated ex with a duplication : ((A:0.1,B:0.2).5:0.3,(C:0.3,(E:0.4,(D_1:0.5,D_2:0.5).D@D:0.5).6:0.05).7:0.3).8:0;
# undated transfer ((B:0.2,A_1:0.1).5:0.3,(C:0.3,(E:0.4,(D:0.5,A_2:0.5).T@D->A:0.5).6:0.05).7:0.3).8:0;

if __name__ == "__main__":
    import sys
    ##############################################
    help = """
Given a file containing reconciled trees in ALE reconciled tree format,
this script writes the trees in recPhyloXML format.

usage : python ALEtoRecPhyloXML.py -g geneFileIn [-o fileOut -s separator]
            -g geneFileIn       : name of the file containing NHX reconciliations
            -o fileOut          : (optional) name of the output file (default is geneFileIn + ".xml" )
            -s separator        : (optional) separator between species and gene name (default: "_")
               """
    #                            (TODO:)
    #                usage : python ALEtoRecPhyloXML.py -g geneFileIn -s speciesFileIn [-o fileOut --include.species]
    #                            (-s speciesFileIn    : (optional) name of the species tree file
    #                            (--include.species   : (optional) whether the species tree should be included in the XML file (using the <spTree> tag)
    
    # infile,ofile = sys.argv[1:]
    
    def parse_umlrec(infile):
        rows = open(infile).read().split("\n")

        stree = rows[2].replace("S:", "").strip()
        spTree = Tree(stree, format=1)
        _s = [idx for idx, _ in enumerate(rows) if "reconciled G-s:" in _][0]
        _e = [
            idx
            for idx, _ in enumerate(rows)
            if "# of\t Duplications\tTransfers\tLosses\tSpeciations" == _
        ][0]
        recgenetrees = rows[_s + 2 : _e]
        return spTree,recgenetrees

    def conv_umlrec2XML(infile,ofile):
        spTree,recgenetrees = parse_umlrec(infile)

        indentLevel = 1
        indentChar = "  "
        OUT = open(ofile, "w")
        OUT.write("<recPhylo>" + "\n")
        spTree = completeTreeNames(spTree, True)
        OUT.write(indentLevel * indentChar + "<spTree>" + "\n")
        indentLevel += 1
        lines = myBasicTreeXMLLines(spTree)
        for xmlline in lines:
            OUT.write(indentLevel * indentChar + xmlline + "\n")
        indentLevel -= 1
        OUT.write(indentLevel * indentChar + "</spTree>" + "\n")
        
        t = recgenetrees[0]
        ALEtree = Tree(t, format=1)
        RT = ALEtreeToReconciledTree(ALEtree, isUndated=True, sepSp="_")
        refineReconciledTreeWithTransferBack(RT)
        ConvertRTtoLossIndepVersion(RT, speciesTree=None, keptChildNameSuffix=".c")
        XMLlines = RT.getTreeRecPhyloXMLLines()
        for xmlline in XMLlines:
            OUT.write(indentLevel * indentChar + xmlline + "\n")
        OUT.write("</recPhylo>" + "\n")
        OUT.close()
        
        
#     # conv_umlrec2XML(infile,ofile)
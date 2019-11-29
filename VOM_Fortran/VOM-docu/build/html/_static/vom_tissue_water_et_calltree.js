function vom_tissue_water_et(node, onCompleteCallback)
{
   var myobj = { label: "vom_mqss", id: "vom_mqss", href: "vom_mqss.html", target:"basefrm" };
   var tmpNode19 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode19.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}

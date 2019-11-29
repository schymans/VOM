function sortcomp(node, onCompleteCallback)
{
   var myobj = { label: "qsort", id: "qsort", href: "qsort.html", target:"basefrm" };
   var tmpNode54 = new YAHOO.widget.TextNode(myobj, node, false);
   tmpNode54.isLeaf = true; 
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}

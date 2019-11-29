function cce(node, onCompleteCallback)
{
   var myobj = { label: "simplex", id: "simplex", href: "simplex.html", target:"basefrm" };
   var tmpNode61 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "sortcomp", id: "sortcomp", href: "sortcomp.html", target:"basefrm" };
   var tmpNode94 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}

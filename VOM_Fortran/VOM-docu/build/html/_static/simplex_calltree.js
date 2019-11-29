function simplex(node, onCompleteCallback)
{
   var myobj = { label: "runmodel", id: "runmodel", href: "runmodel.html", target:"basefrm" };
   var tmpNode74 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}

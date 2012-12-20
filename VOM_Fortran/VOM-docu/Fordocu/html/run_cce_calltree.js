function run_cce(node, onCompleteCallback)
{
   var myobj = { label: "cce", id: "cce", href: "cce.html", target:"basefrm" };
   var tmpNode72 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}

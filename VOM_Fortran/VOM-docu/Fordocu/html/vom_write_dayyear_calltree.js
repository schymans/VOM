function vom_write_dayyear(node, onCompleteCallback)
{
   var myobj = { label: "vom_add_yearly", id: "vom_add_yearly", href: "vom_add_yearly.html", target:"basefrm" };
   var tmpNode45 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}

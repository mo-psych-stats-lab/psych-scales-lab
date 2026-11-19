/* ===== Privacy agreement (shared by both labs) =====
   Signed-in students see this once per account and must accept it. After acceptance their results from
   «اعرف نفسك» and «مكتبة الاختبارات» are saved to the lecturer's admin panel (written consent is also
   collected on paper). Nothing is published or shared with anyone else. ===== */
(function () {
  'use strict';
  var KEY = 'privacy-v1';
  var get = function () { try { return JSON.parse(localStorage.getItem(KEY) || '{}') || {}; } catch (e) { return {}; } };
  var fsOK = function () { try { return !!(window.ACCT && ACCT.fs && ACCT.fs.ready()); } catch (e) { return false; } };
  var uid = function () { return fsOK() ? ACCT.fs.uid() : null; };
  function ok() { var u = uid(); return !!(u && get()[u]); }
  var TEXT = '<h3 style="margin:0">سياسة الخصوصية والموافقة</h3>' +
    '<p style="margin:0">قبل ما تكمل، اقرأ إزاي بنتعامل مع بياناتك:</p><ul style="margin:0;line-height:1.9;padding-inline-start:20px">' +
    '<li><b>اللي بيتحفظ:</b> اسمك وإيميلك، وتقدمك في الدروس والتدريب، ونتايج «اعرف نفسك» و«مكتبة الاختبارات».</li>' +
    '<li><b>مين يشوفها:</b> مسؤول الموقع (المحاضر) بس، علشان المتابعة الأكاديمية وتوجيهك لو احتجت مساعدة. ومفيش طالب يقدر يشوف بيانات طالب تاني.</li>' +
    '<li><b>مش بتتنشر:</b> بياناتك مش بتتنشر ولا بتتباع ولا بتتشارك مع أي جهة تانية. ولو اتستخدمت في بحث بتكون مجمّعة ومن غير أسماء.</li>' +
    '<li><b>النتايج النفسية فرز تعليمي مش تشخيص،</b> ولو نتيجتك محتاجة اهتمام هتلاقي أرقام خطوط الدعم في صفحة النتائج.</li>' +
    '<li><b>حقك:</b> تقدر تطلب حذف بياناتك في أي وقت، أو تحذف حسابك بنفسك من «حسابي».</li></ul>' +
    '<label style="display:flex;gap:8px;align-items:flex-start"><input type="checkbox" id="pvChk"> <span>قريت وموافق على سياسة الخصوصية دي.</span></label>' +
    '<div style="display:flex;gap:8px;flex-wrap:wrap"><button type="button" id="pvYes" style="border:0;border-radius:10px;padding:8px 18px;font:inherit;font-weight:700;background:var(--primary,var(--teal,#4f46e5));color:#fff;cursor:pointer">موافق</button><button type="button" id="pvNo" style="border:1px solid currentColor;border-radius:10px;padding:8px 18px;font:inherit;background:transparent;color:inherit;cursor:pointer">مش دلوقتي</button></div><p id="pvMsg" style="margin:0;font-size:13px"></p>';
  var shown = false;
  function show() {
    if (shown || document.getElementById('pvOv')) return; shown = true;
    var ov = document.createElement('div'); ov.id = 'pvOv'; ov.setAttribute('role', 'dialog'); ov.setAttribute('aria-modal', 'true'); ov.setAttribute('aria-labelledby', 'pvT');
    ov.style.cssText = 'position:fixed;inset:0;z-index:9999;background:rgba(0,0,0,.45);display:grid;place-items:center;padding:16px';
    ov.innerHTML = '<div dir="rtl" style="max-width:560px;width:100%;max-height:90vh;overflow:auto;background:var(--card,var(--surface,#fff));color:var(--ink,inherit);border-radius:16px;padding:18px;display:flex;flex-direction:column;gap:10px;line-height:1.8;box-shadow:0 10px 40px rgba(0,0,0,.25)">' + TEXT + '</div>';
    document.body.appendChild(ov);
    var chk = ov.querySelector('#pvChk'); chk.focus();
    ov.querySelector('#pvYes').onclick = function () {
      if (!chk.checked) { ov.querySelector('#pvMsg').textContent = 'علّم على «قريت وموافق» الأول.'; return; }
      var u = uid(), m = get(); if (u) { m[u] = Date.now(); try { localStorage.setItem(KEY, JSON.stringify(m)); } catch (e) { } }
      ov.remove(); window.dispatchEvent(new Event('privacy-accepted'));
    };
    ov.querySelector('#pvNo').onclick = function () { ov.remove(); };
  }
  var last = null, fired = {};
  (function poll(n) {
    var u = uid();
    if (u !== last) { last = u; shown = false; }
    if (u && !location.hash.startsWith('#admin')) { if (!get()[u]) show(); else if (!fired[u]) { fired[u] = 1; window.dispatchEvent(new Event('privacy-accepted')); } }
    setTimeout(function () { poll(n + 1); }, n < 20 ? 1500 : 5000);
  })(0);
  document.addEventListener('click', function (e) { var t = e.target.closest && e.target.closest('[data-pvshow]'); if (t) { e.preventDefault(); shown = false; show(); } });
  window.PRIVACY = { ok: ok, show: function () { shown = false; show(); } };
})();

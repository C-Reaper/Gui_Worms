#include "/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
#include "/home/codeleaded/System/Static/Library/TransformedView.h"


#define STEP_SIZE       0.1f
#define POSITION_SIZE   0.5f

float angle;
float strength;
Vec2 position;
Vec2 force;
Vec2 target;

Vec2* focused;
TransformedView tv;

// f(t) = src + d(a) * s * t + 0.5f * g * t * t
// dst = f(t)
// 0 = src + d(a) * s * t + 0.5f * g * t * t - dst
// I:   0 = src.x - dst.x + cos(a) * s * t + 0.5f * g.x * t * t
// II:  0 = src.y - dst.y + sin(a) * s * t + 0.5f * g.y * t * t
//
// t1 = (-cos(a) * s + sqrt((cos(a) * s)^2 - 4 * g.x * (src.x - dst.x))) / 2 * g.x
// t2 = (-cos(a) * s - sqrt((cos(a) * s)^2 - 4 * g.x * (src.x - dst.x))) / 2 * g.x
// t1,2 in II:
// t1: 0 = src.y - dst.y + sin(a) * s * t1 + 0.5f * g.y * t1 * t1
// t2: 0 = src.y - dst.y + sin(a) * s * t2 + 0.5f * g.y * t2 * t2
// => Newton: xn+1 = xn - f(xn) / f'(xn)
float Angle_FuncT1(Vec2 src,Vec2 dst,Vec2 g,float s,float a){
    const float c_1 = src.y - dst.y;
    const float b_1 = sinf(a) * s;
    const float a_1 = 0.5f * g.y;

    const float abc_c = src.x - dst.x;
    const float abc_b = cosf(a) * s;
    const float abc_a = g.x;
    const float ds = abc_b * abc_b - 4.0f * abc_a * abc_c;
    if(ds < 0.0f) return NAN;
    
    const float t = (-abc_b + sqrtf(ds)) / (2.0f * abc_a);
    
    return c_1 + b_1 * t + a_1 * t * t;
}
float Angle_FuncPrimeT1(Vec2 src,Vec2 dst,Vec2 g,float s,float a){
    const float a_1 = g.y;

    const float abc_c = src.x - dst.x;
    const float abc_b = cosf(a) * s;
    const float abc_a = g.x;

    const float bp  = -s * sinf(a);

    const float ds = abc_b * abc_b - 4.0f * abc_a * abc_c;
    if(ds < 0.0f) return NAN;

    const float dp  = 2.0f * abc_b * bp;

    const float t = (-abc_b + sqrtf(ds)) / (2.0f * abc_a);
    const float tp = (-bp + dp / (2.0f * sqrtf(ds))) / (2.0f * abc_a);
    
    return s * (cosf(a) * t + sinf(a) * tp) + a_1 * t * tp;
}
float Angle_FuncT2(Vec2 src,Vec2 dst,Vec2 g,float s,float a){
    const float c_1 = src.y - dst.y;
    const float b_1 = sinf(a) * s;
    const float a_1 = 0.5f * g.y;

    const float abc_c = src.x - dst.x;
    const float abc_b = cosf(a) * s;
    const float abc_a = g.x;
    const float ds = abc_b * abc_b - 4.0f * abc_a * abc_c;
    if(ds < 0.0f) return NAN;
    
    const float t = (-abc_b - sqrtf(ds)) / (2.0f * abc_a);
    
    return c_1 + b_1 * t + a_1 * t * t;
}
float Angle_FuncPrimeT2(Vec2 src,Vec2 dst,Vec2 g,float s,float a){
    const float a_1 = g.y;

    const float abc_c = src.x - dst.x;
    const float abc_b = cosf(a) * s;
    const float abc_a = g.x;

    const float bp  = -s * sinf(a);

    const float ds = abc_b * abc_b - 4.0f * abc_a * abc_c;
    if(ds < 0.0f) return NAN;

    const float dp  = 2.0f * abc_b * bp;

    const float t = (-abc_b - sqrtf(ds)) / (2.0f * abc_a);
    const float tp = (-bp - dp / (2.0f * sqrtf(ds))) / (2.0f * abc_a);
    
    return s * (cosf(a) * t + sinf(a) * tp) + a_1 * t * tp;
}

int Angle_Find_Newton(Vec2 src,Vec2 dst,Vec2 g,float s,float* a0,float (*f)(Vec2,Vec2,Vec2,float,float),float (*fp)(Vec2,Vec2,Vec2,float,float)){
    float anp1 = *a0;
    for(int i = 0;i<100;i++){
        float out_f = f(src,dst,g,s,anp1);
        float out_fp = fp(src,dst,g,s,anp1);
        if (fabsf(out_fp) < 1e-5f) return 0;
        
        const float newanp1 = anp1 - out_f / out_fp;
        //printf("%d: %f\t(%f,%f)\t-> %f\n",i,anp1,out_f,out_fp,newanp1);
        //printf("%d: %f -> %f\n",i,anp1,newanp1);
        if(isnan(out_f) || isnan(out_fp)) return 0;
        if (fabsf(newanp1 - anp1) < 1e-6f) break;
        
        anp1 = newanp1;
    }
    *a0 = anp1;
    return 1;
}
int Angle_Calc(Vec2 src,Vec2 dst,Vec2 g,float s,float* a1,float* a2){
    //const int found1 = Angle_Find_Newton(src,dst,g,s,a1,Angle_FuncT1,Angle_FuncPrimeT1);
    //const int found2 = Angle_Find_Newton(src,dst,g,s,a2,Angle_FuncT2,Angle_FuncPrimeT2);
    //return found1 + found2;

    float fa = Vec2_AngleOf(g) - F32_PI05;

    Vec2 nrd = Vec2_Sub(dst,src);
    Vec2 d = Vec2_Mulf(Vec2_OfAngle(Vec2_AngleOf(nrd) - fa),Vec2_Mag(nrd));
    
    float a = 0.5f * Vec2_Mag(g) * (d.x / s) * (d.x / s);
    float b = d.x;
    float c = 0.5f * Vec2_Mag(g) * (d.x / s) * (d.x / s) - d.y;

    float dis = b * b - 4.0f * a * c;

    float dir = (d.x >= 0.0f ? 0.0f : F32_PI) + fa;
    *a1 = dir + atanf((-b + sqrtf(dis)) / (2.0f * a));
    *a2 = dir + atanf((-b - sqrtf(dis)) / (2.0f * a));

    if(dis < 0.0f)    return 0;
    if(dis == 0.0f)   return 1;
    return 2;
}

void Trail_Render(unsigned int* Target,int Target_Width,int Target_Height,Vec2 p,float a,float s,unsigned int c){
    const Vec2 startv = Vec2_Mulf(Vec2_OfAngle(a),s);

    for(float i = 0.0f;i<10.0f;i+=STEP_SIZE){
        const Vec2 dir = Vec2_Add(startv,Vec2_Mulf(force,0.5f * i));
        const Vec2 pos = Vec2_Add(p,Vec2_Mulf(dir,i));
        const Vec2 out = TransformedView_WorldScreenPos(&tv,pos);
        Point_RenderX(WINDOW_STD_ARGS,out.x,out.y,c);
    }
}

void Setup(AlxWindow* w){
    angle = -F32_PI025;
    strength = 1.0f;
    position = (Vec2){ 0.0f,0.0f };
    force = (Vec2){ 1.0f,1.0f };
    target = (Vec2){ 20.0f,1.0f };

    focused = NULL;
    tv = TransformedView_New((Vec2){ GetWidth(),GetHeight() });
    TransformedView_Zoom(&tv,(Vec2){ 0.05f,0.05f });
}
void Update(AlxWindow* w){
    TransformedView_Output(&tv,(Vec2){ GetWidth(),GetHeight() });
    TransformedView_HandlePanZoom(&tv,w->Strokes,GetMouse());

    if(Stroke(ALX_KEY_LEFT).DOWN){
        angle -= F32_PI * w->ElapsedTime;
    }
    if(Stroke(ALX_KEY_RIGHT).DOWN){
        angle += F32_PI * w->ElapsedTime;
    }
    if(Stroke(ALX_KEY_R).DOWN){
        strength -= 1.0f * w->ElapsedTime;
    }
    if(Stroke(ALX_KEY_F).DOWN){
        strength += 1.0f * w->ElapsedTime;
    }

    Circle p_r = Circle_New(position,POSITION_SIZE);
    Circle p_t = Circle_New(target  ,POSITION_SIZE);
    const Vec2 m_w = TransformedView_ScreenWorldPos(&tv,GetMouse()); 
    if(Stroke(ALX_MOUSE_L).PRESSED){
        if(Circle_Point(&p_r,m_w))
            focused = &position;
        if(Circle_Point(&p_t,m_w))
            focused = &target;
    }
    if(Stroke(ALX_MOUSE_L).RELEASED){
        focused = NULL;
    }
    if(focused) *focused = m_w;

    Clear(BLACK);

    if(Stroke(ALX_KEY_T).DOWN){
        float a1;// = Vec2_AngleOf(Vec2_Sub(Vec2_Sub(target,position),force));
        float a2;// = Vec2_AngleOf(Vec2_Sub(Vec2_Sub(target,position),force));
        const int count = Angle_Calc(position,target,force,strength,&a1,&a2);

        if(count > 0) Trail_Render(WINDOW_STD_ARGS,position,a1,strength,YELLOW);
        if(count > 1) Trail_Render(WINDOW_STD_ARGS,position,a2,strength,LIGHT_YELLOW);
    }
    
    Circle_RenderX(WINDOW_STD_ARGS,TransformedView_WorldScreenPos(&tv,p_r.p),TransformedView_WorldScreenLX(&tv,p_r.r),BLUE);
    Circle_RenderX(WINDOW_STD_ARGS,TransformedView_WorldScreenPos(&tv,p_t.p),TransformedView_WorldScreenLX(&tv,p_t.r),GREEN);

    Trail_Render(WINDOW_STD_ARGS,position,angle,strength,RED);
}
void Delete(AlxWindow* w){

}

int main() {
    if(Create("Worms",1800,1300,1,1,Setup,Update,Delete)){
        Start();
    }
    return 0;
}
